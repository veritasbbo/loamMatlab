function [pc_sharpXyzIdRt, pc_lessSharpXyzIdRt, pc_flatXyzIdRt, pc_lessFlatXyzIdRt, pcAlignedXyzIdRt, time] = scanRegistration(pcRawXyz, time, N_SCANS, MINIMUM_RANGE)
    
    t_pre = tic;
    %% 
    pcRawXyz = removeClosedPointCloud(pcRawXyz, MINIMUM_RANGE);
    
    %% 스캔 라인에 벗어난 측정치 제거
    angle = atan(pcRawXyz(:,3)./sqrt(pcRawXyz(:,1).^2 + pcRawXyz(:,2).^2))*180/pi;
    if(N_SCANS == 16)
        scanId = floor((angle + 15)./2 + 0.5);
        outlierInd = or(scanId > (N_SCANS - 1), scanId < 0);
    elseif(N_SCANS == 32)
        scanId = floor((angle + 92/3)*3/4);
        outlierInd = or(scanId > (N_SCANS - 1), scanId < 0);
    elseif(N_SCANS == 64)
        scanId = nan(length(angle),1);
        scanId(angle>=-8.83) = floor((2-angle(angle>=-8.83))*3 + 0.5);
        scanId(angle<-8.83) = N_SCANS/2 + floor((-8.83 - angle(angle<-8.83))*2 + 0.5);
        outlierInd = or(or(or(scanId > 50, scanId<0),angle<-24.33),angle>2);
    end
    cloudSize = sum(~outlierInd);
    scanId(outlierInd) = [];
    scanId = scanId + 1;    % matlab 1부터 시작 반영
    pcRawXyz(outlierInd) = [];
    
    %% 포인트별 상대 시간
    ori = -atan2(pcRawXyz(:,2), pcRawXyz(:,1));
    startOri = ori(1);
    ori = ori - startOri;
    endOri = ori(end) + 2*pi;
    if(endOri > 3*pi)
        endOri = endOri - 2*pi;
    end
    if(endOri < pi)
        endOri = endOri + 2*pi;
    end
    halfPassed = false;
    for i=1:cloudSize
        if(~halfPassed)
            if(ori(i)<- pi/2)
                ori(i) = ori(i) + 2*pi;
            elseif(ori(i)> pi*3/2)
                ori(i) = ori(i) - 2*pi;
            end
            if(ori(i)>pi)
                halfPassed = true;
                ori(i+1:end) = ori(i+1:end) + 2*pi;
            end
        else
            if(ori(i)<endOri - pi*3/2)
                ori(i) = ori(i) + 2*pi;
            elseif(ori(i)>endOri + pi/2)
                ori(i) = ori(i) - 2*pi;
            end
        end
    end
    relTimeRto = ori/ori(end);
    relTime = relTimeRto;
    %% 스캔ID 순서대로 포인트 정리
    scanIdAlignedIndTemp = [1:cloudSize]';
    scanIdAlignedInd = nan(cloudSize,1);
    scanStartInd = zeros(N_SCANS,1);
    scanEndInd = zeros(N_SCANS,1);
    regLen = 0;
    for scanIdx = 1:N_SCANS
        scanStartInd(scanIdx) = regLen + 6;
        pcId4curScanId = scanIdAlignedIndTemp(scanId == scanIdx);
        regLenNew = regLen + length(pcId4curScanId);
        scanEndInd(scanIdx) = regLenNew - 5;
        scanIdAlignedInd(regLen + 1:regLenNew) = pcId4curScanId;
        regLen = regLenNew;
    end
    pcAlignedXyzIdRt = [pcRawXyz(scanIdAlignedInd,:),scanId(scanIdAlignedInd),relTime(scanIdAlignedInd,:)];
%     pcAlignedRng = sqrt(sum(pcAlignedXyz.^2,2));
    time.t_prepare = time.t_prepare + toc(t_pre);
    %% cuvature 정리
    diffXyz = zeros(cloudSize - 10,3);
    diffXyz(:,1) = pcAlignedXyzIdRt(1:end - 10,1) + pcAlignedXyzIdRt(2:end - 9,1) + pcAlignedXyzIdRt(3:end - 8,1) + pcAlignedXyzIdRt(4:end - 7,1) + pcAlignedXyzIdRt(5:end - 6,1) ...
            - pcAlignedXyzIdRt(6:end - 5,1)*10 ...
            + pcAlignedXyzIdRt(7:end - 4,1) + pcAlignedXyzIdRt(8:end - 3,1) + pcAlignedXyzIdRt(9:end - 2,1) + pcAlignedXyzIdRt(10:end - 1,1) + pcAlignedXyzIdRt(11:end,1);
    diffXyz(:,2) = pcAlignedXyzIdRt(1:end - 10,2) + pcAlignedXyzIdRt(2:end - 9,2) + pcAlignedXyzIdRt(3:end - 8,2) + pcAlignedXyzIdRt(4:end - 7,2) + pcAlignedXyzIdRt(5:end - 6,2) ...
            - pcAlignedXyzIdRt(6:end - 5,2)*10 ...
            + pcAlignedXyzIdRt(7:end - 4,2) + pcAlignedXyzIdRt(8:end - 3,2) + pcAlignedXyzIdRt(9:end - 2,2) + pcAlignedXyzIdRt(10:end - 1,2) + pcAlignedXyzIdRt(11:end,2);
    diffXyz(:,3) = pcAlignedXyzIdRt(1:end - 10,3) + pcAlignedXyzIdRt(2:end - 9,3) + pcAlignedXyzIdRt(3:end - 8,3) + pcAlignedXyzIdRt(4:end - 7,3) + pcAlignedXyzIdRt(5:end - 6,3) ...
            - pcAlignedXyzIdRt(6:end - 5,3)*10 ...
            + pcAlignedXyzIdRt(7:end - 4,3) + pcAlignedXyzIdRt(8:end - 3,3) + pcAlignedXyzIdRt(9:end - 2,3) + pcAlignedXyzIdRt(10:end - 1,3) + pcAlignedXyzIdRt(11:end,3);
    cloudCurvature = [zeros(5,1);sum(diffXyz.^2,2);zeros(5,1)];
    cloudLabel = [nan(5,1);zeros(cloudSize - 10,1);nan(5,1)];
    %%
    t_pts = tic;
    pc_sharpXyzIdRt = nan(cloudSize,5);
    pc_lessSharpXyzIdRt = nan(cloudSize,5);
    pc_flatXyzIdRt = nan(cloudSize,5);
    pc_lessFlatXyzIdRt = nan(cloudSize,5);
    
    cloudNeighborPicked = false(cloudSize,1);
    cuvatureSortedInd = [1:cloudSize]';
    
    sharpCloudSize = 0;
    lessSharpCloudSize = 0;
    flatCloudSize = 0;
    lessFlatCloudSize = 0;
    % 6개 구역 분할
    for scanIdx = 1:N_SCANS
        if(scanEndInd(scanIdx) - scanStartInd(scanIdx) < 6)
           continue;
        end
        for subsecIdx = 0:5
            sp = floor(scanStartInd(scanIdx) + (scanEndInd(scanIdx) - scanStartInd(scanIdx))*subsecIdx/6);
            ep = floor(scanStartInd(scanIdx) + (scanEndInd(scanIdx) - scanStartInd(scanIdx))*(subsecIdx + 1)/6 - 1);
            [~, sortedLocalInd] = sort(cloudCurvature(sp:ep));
            cuvatureSortedInd(sp:ep) = cuvatureSortedInd(sp - 1 + sortedLocalInd,:);

            % pc_sharpXyz 추출
            largestPickedNum = 0;
            for i = ep:-1:sp
        %                diffXyz0 = [1,0;0,-1]*diff(pcAlignedXyz(cuvatureSortedInd(i) - 1:cuvatureSortedInd(i) + 1,:),1,1);
        %                magDiffXyz = sqrt(sum(diffXyz0.^2,2));
        %                obscuredTest = any(pcAlignedXyz(cuvatureSortedInd(i),:)*diffXyz0'./magDiffXyz'/pcAlignedRng(cuvatureSortedInd(i))>0.95);
        %                if(cloudCurvature(cuvatureSortedInd(i)) > 0.1 && ~cloudNeighborPicked(cuvatureSortedInd(i))&&~obscuredTest)
                if(cloudCurvature(cuvatureSortedInd(i)) > 0.1 && ~cloudNeighborPicked(cuvatureSortedInd(i)))
                largestPickedNum = largestPickedNum + 1;
                if (largestPickedNum <= 2)
                    sharpCloudSize = sharpCloudSize + 1;
                    lessSharpCloudSize = lessSharpCloudSize + 1;
                    pc_sharpXyzIdRt(sharpCloudSize,:) = pcAlignedXyzIdRt(cuvatureSortedInd(i),:);
                    pc_lessSharpXyzIdRt(lessSharpCloudSize,:) = pcAlignedXyzIdRt(cuvatureSortedInd(i),:);
                    cloudLabel(cuvatureSortedInd(i)) = 2;
                elseif (largestPickedNum <= 20)
                    lessSharpCloudSize = lessSharpCloudSize + 1;
                    pc_lessSharpXyzIdRt(lessSharpCloudSize,:) = pcAlignedXyzIdRt(cuvatureSortedInd(i),:);
                    cloudLabel(cuvatureSortedInd(i)) = 1;
                else
                    break;
                end
                end
                cloudNeighborPicked(cuvatureSortedInd(i)) = true;

                diffXyz = pcAlignedXyzIdRt(cuvatureSortedInd(i) + 1:cuvatureSortedInd(i) + 5,:) - pcAlignedXyzIdRt(cuvatureSortedInd(i):cuvatureSortedInd(i) + 4,:);
                test = find(sum(diffXyz.^2,2) > 0.05);
                if(~isempty(test))
                    neighborPicked = test(1);
                    cloudNeighborPicked(cuvatureSortedInd(i) + 1:cuvatureSortedInd(i) + neighborPicked - 1,:) = true;
                end
                diffXyz = pcAlignedXyzIdRt(cuvatureSortedInd(i) - 1:-1:cuvatureSortedInd(i) - 5,:) - pcAlignedXyzIdRt(cuvatureSortedInd(i):-1:cuvatureSortedInd(i) - 4,:);
                test = find(sum(diffXyz.^2,2) > 0.05);
                if(~isempty(test))
                    neighborPicked = test(1);
                    cloudNeighborPicked(cuvatureSortedInd(i) - 1:-1:cuvatureSortedInd(i) - neighborPicked + 1,:) = true;
                end
           end

           % pc_flatXyz 추출
           smallestPickedNum = 0;
           for i = sp:1:ep
        %                diffXyz0 = diff(pcAlignedXyz(cuvatureSortedInd(i) - 1:cuvatureSortedInd(i) + 1,:),1,1);
        %                magDiffXyz = sqrt(sum(diffXyz0.^2,2));
        %                parTest = any(abs(pcAlignedXyz(cuvatureSortedInd(i),:)*diffXyz0'./magDiffXyz'/pcAlignedRng(cuvatureSortedInd(i)))>0.95);
        %                if(cloudCurvature(cuvatureSortedInd(i)) < 0.1 && ~cloudNeighborPicked(cuvatureSortedInd(i)) && ~parTest)
               if(cloudCurvature(cuvatureSortedInd(i)) < 0.1 && ~cloudNeighborPicked(cuvatureSortedInd(i)))
                    smallestPickedNum = smallestPickedNum + 1;
                    if (smallestPickedNum <= 4)
                        flatCloudSize = flatCloudSize + 1;
                        pc_flatXyzIdRt(flatCloudSize,:) = pcAlignedXyzIdRt(cuvatureSortedInd(i),:);
                        cloudLabel(cuvatureSortedInd(i)) = -1;
                    else
                        break;
                    end
               end
               cloudNeighborPicked(cuvatureSortedInd(i)) = true;

               diffXyz = pcAlignedXyzIdRt(cuvatureSortedInd(i) + 1:cuvatureSortedInd(i) + 5,:) - pcAlignedXyzIdRt(cuvatureSortedInd(i):cuvatureSortedInd(i) + 4,:);
               test = find(sum(diffXyz.^2,2) > 0.05);
               if(~isempty(test))
                   neighborPicked = test(1);
                   cloudNeighborPicked(cuvatureSortedInd(i) + 1:cuvatureSortedInd(i) + neighborPicked - 1,:) = true;
               end
               diffXyz = pcAlignedXyzIdRt(cuvatureSortedInd(i) - 1:-1:cuvatureSortedInd(i) - 5,:) - pcAlignedXyzIdRt(cuvatureSortedInd(i):-1:cuvatureSortedInd(i) - 4,:);
               test = find(sum(diffXyz.^2,2) > 0.05);
               if(~isempty(test))
                   neighborPicked = test(1);
                   cloudNeighborPicked(cuvatureSortedInd(i) - 1:-1:cuvatureSortedInd(i) - neighborPicked + 1,:) = true;
               end
           end

        end
        pcScanAlignedXyzIdRt = pcAlignedXyzIdRt(scanStartInd(scanIdx):scanEndInd(scanIdx),:);
        surfPointsLessFlatScanXyzIdRt = pcScanAlignedXyzIdRt(cloudLabel(scanStartInd(scanIdx):scanEndInd(scanIdx))<=0,:);
        surfPointsLessFlatScanCloud = pointCloud(surfPointsLessFlatScanXyzIdRt(:,1:3),'Intensity',surfPointsLessFlatScanXyzIdRt(:,5));
        pcScan_lessFlatCloud = pcdownsample(surfPointsLessFlatScanCloud, 'gridAverage', 0.2);
        pcScan_lessFlatXyzIdRt = [pcScan_lessFlatCloud.Location,scanIdx*ones(pcScan_lessFlatCloud.Count,1), pcScan_lessFlatCloud.Intensity];
        pc_lessFlatXyzIdRt(lessFlatCloudSize + 1: lessFlatCloudSize + pcScan_lessFlatCloud.Count,:) = pcScan_lessFlatXyzIdRt;
        lessFlatCloudSize = lessFlatCloudSize + pcScan_lessFlatCloud.Count;
    end

    pc_sharpXyzIdRt(sharpCloudSize + 1:end,:) = [];
    pc_lessSharpXyzIdRt(lessSharpCloudSize + 1:end,:) = [];
    pc_flatXyzIdRt(flatCloudSize + 1:end,:) = [];
    pc_lessFlatXyzIdRt(lessFlatCloudSize + 1:end,:) = [];
    time.t_pts = time.t_pts + toc(t_pts);
end