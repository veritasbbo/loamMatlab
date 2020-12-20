function [pc_cornerLastXyzIdRt, pc_surfLastXyzIdRt, pcAlignedXyzIdRt, qt_w_curr, qt_last_curr, timeVar, systemInited] = laserOdometry(pc_sharpXyzIdRt, pc_lessSharpXyzIdRt, pc_flatXyzIdRt, pc_lessFlatXyzIdRt, pc_cornerLastXyzIdRt, pc_surfLastXyzIdRt, pcAlignedXyzIdRt, qt_w_curr, qt_last_curr, timeVar, systemInited, SCAN_LINE, DISTANCE_SQ_THRESHOLD, NEARBY_SCAN, HUBER_THRES)
    
    q_w_curr = qt_w_curr(1:4);
    t_w_curr = qt_w_curr(5:7);
    q_last_curr = qt_last_curr(1:4);
    t_last_curr = qt_last_curr(5:7);
    DISTORTION = false;
       
    if(~systemInited)
        systemInited = true;
    else
        t_opt = tic;
%         startEndInd4cornerfLast = zeros(SCAN_LINE,2);
%         startEndInd4surfLast = zeros(SCAN_LINE,2);
%         for scanId = 1:SCAN_LINE
%             indTmp = find(pc_cornerLastXyzIdRt(:,4)==scanId);
%             startEndInd4cornerfLast(scanId,:) = [min(indTmp),max(indTmp)];
%             indTmp = find(pc_surfLastXyzIdRt(:,4)==scanId);
%             startEndInd4surfLast(scanId,:) = [min(indTmp),max(indTmp)];
%         end
        for opti_counter = 0:1
            %%
            t_data = tic;
            %% find correspondence for corner features
            pointSelXyzIdRt = transformToStart(pc_sharpXyzIdRt, q_last_curr, t_last_curr, DISTORTION);
            [closestPointInd, D] = knnsearch(pc_cornerLastXyzIdRt(:,1:3),pointSelXyzIdRt(:,1:3),'K',1,'NSMethod','kdtree');
            innerTest = D < DISTANCE_SQ_THRESHOLD;
            innerQueryPoints = pc_sharpXyzIdRt(innerTest,:);
            innerDistclosestPointInd = closestPointInd(innerTest);
            innerDistPoints = pc_cornerLastXyzIdRt(innerDistclosestPointInd, :);
            numOfInnerDistPoints = length(innerDistPoints(:,1));
            minPointArr = zeros(numOfInnerDistPoints,3);
            minPointSqDis2 = zeros(numOfInnerDistPoints,1);
            
            for corrId = 1:numOfInnerDistPoints
                closestPointScanID = pc_cornerLastXyzIdRt(innerDistclosestPointInd(corrId),4);
                nearScanPointsSelector = and(pc_cornerLastXyzIdRt(:,4) ~= closestPointScanID, abs(pc_cornerLastXyzIdRt(:,4) - closestPointScanID) < NEARBY_SCAN);
                nearScanPoints = pc_cornerLastXyzIdRt(nearScanPointsSelector,:);
                [minPointInd2, minPointSqDis2(corrId)] = knnsearch(nearScanPoints(:,1:3),innerDistPoints(corrId,1:3),'K',1,'NSMethod','exhaustive');
                minPointArr(corrId,:) = nearScanPoints(minPointInd2,1:3);
            end
            hasNear2Selector = minPointSqDis2 < DISTANCE_SQ_THRESHOLD;
            corner_curr_point = innerQueryPoints(hasNear2Selector,:);
            corner_last_point_a = pc_cornerLastXyzIdRt(innerDistclosestPointInd(hasNear2Selector),1:3);
            corner_last_point_b = minPointArr(hasNear2Selector,:);
            numCornerCorresponds = sum(hasNear2Selector);
            %% find correspondence for plane features
            pointSelXyzIdRt = transformToStart(pc_flatXyzIdRt, q_last_curr, t_last_curr, DISTORTION);
            [closestPointInd, D] = knnsearch(pc_surfLastXyzIdRt(:,1:3),pointSelXyzIdRt(:,1:3),'K',1,'NSMethod','kdtree');
            innerTest = D < DISTANCE_SQ_THRESHOLD;
            innerQueryPoints = pc_flatXyzIdRt(innerTest,:);
            innerDistclosestPointInd = closestPointInd(innerTest);
            innerDistPoints = pc_surfLastXyzIdRt(innerDistclosestPointInd, :);
            numOfInnerDistPoints = length(innerDistPoints(:,1));
            minPointSqDis2 = zeros(numOfInnerDistPoints,1);
            minPointSqDis3 = zeros(numOfInnerDistPoints,1);
            minPoint2Arr = zeros(numOfInnerDistPoints,3);
            minPoint3Arr = zeros(numOfInnerDistPoints,3);
            for corrId = 1:numOfInnerDistPoints
                closestPointScanID = pc_surfLastXyzIdRt(innerDistclosestPointInd(corrId),4);
                pc_surfLastXyzIdRtExceptClosestPoint = pc_surfLastXyzIdRt([1:innerDistclosestPointInd(corrId)-1,innerDistclosestPointInd(corrId)+1:end],:);
                sameScanPoints = pc_surfLastXyzIdRtExceptClosestPoint(pc_surfLastXyzIdRtExceptClosestPoint(:,4) == closestPointScanID,:);
                [minPointInd2, minPointSqDis2(corrId)] = knnsearch(sameScanPoints(:,1:3),innerDistPoints(corrId,1:3),'K',1,'NSMethod','exhaustive');
                minPoint2Arr(corrId,:) = sameScanPoints(minPointInd2,1:3);
                nearScanPointsSelector = and(pc_surfLastXyzIdRt(:,4) ~= closestPointScanID, abs(pc_surfLastXyzIdRt(:,4) - closestPointScanID) < NEARBY_SCAN);
                nearScanPoints = pc_surfLastXyzIdRt(nearScanPointsSelector,:);
                [minPointInd3, minPointSqDis3(corrId)] = knnsearch(nearScanPoints(:,1:3),innerDistPoints(corrId,1:3),'K',1,'NSMethod','exhaustive');
                minPoint3Arr(corrId,:) = nearScanPoints(minPointInd3,1:3);
            end
            hasNear23Selector = and(minPointSqDis2 < DISTANCE_SQ_THRESHOLD, minPointSqDis3 < DISTANCE_SQ_THRESHOLD);
            
            surf_curr_point = innerQueryPoints(hasNear23Selector,:);
            surf_last_point_a = pc_surfLastXyzIdRt(innerDistclosestPointInd(hasNear23Selector),1:3);
            surf_last_point_b = minPoint2Arr(hasNear23Selector,:);
            surf_last_point_c = minPoint3Arr(hasNear23Selector,:);
            numSurfCorresponds = sum(hasNear23Selector);
            timeVar.t_data = timeVar.t_data + toc(t_data);
            %%
            if(numCornerCorresponds + numSurfCorresponds < 10)
                warning('less correspondence!');
            end
            t_solver = tic;
            fun = @(x)residualFunc(x, corner_curr_point, corner_last_point_a, corner_last_point_b, surf_curr_point, surf_last_point_a, surf_last_point_b,surf_last_point_c, DISTORTION, HUBER_THRES);
            options.Algorithm = 'levenberg-marquardt';
            options.Display = 'final-detailed';
            options.Display = 'off';
            options.MaxIterations = 4;
            qt_last_curr = [q_last_curr;t_last_curr];
            qt_last_curr = lsqnonlin(fun,qt_last_curr,[],[],options);
            q_last_curr = qt_last_curr(1:4)/norm(qt_last_curr(1:4));
            t_last_curr = qt_last_curr(5:7);
            timeVar.t_solver = timeVar.t_solver + toc(t_solver);
            timeVar.t_opt = timeVar.t_opt + toc(t_opt);
        end
    end
    
    t_w_curr = t_w_curr + rotateframe(quaternion(q_w_curr'),t_last_curr')';
    q_w_curr = compact(normalize(quaternion(q_w_curr')*quaternion(q_last_curr')))';
    
    qt_w_curr = [q_w_curr;t_w_curr];
    qt_last_curr = [q_last_curr;t_last_curr];
    % transform corner features and plane features to the scan end point
    if (false)
        pc_lessSharpXyzIdRt = transformToEnd(pc_lessSharpXyzIdRt);
        pc_lessFlatXyzIdRt = transformToEnd(pc_lessFlatXyzIdRt);
        pcAlignedXyzIdRt = transformToEnd(pcAlignedXyzIdRt);
    end
    pc_cornerLastXyzIdRt = pc_lessSharpXyzIdRt;
    pc_surfLastXyzIdRt = pc_lessFlatXyzIdRt;
end