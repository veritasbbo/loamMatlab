%% Parameters
SCAN_LINE = 16;
SCAN_PERIOD = 0.1;
DISTANCE_SQ_THRESHOLD = 25;
NEARBY_SCAN = 2.5;
MAPPING_SKIP_FRAME = 1;
MINIMUM_RANGE = 0.3;
MAPPING_LINE_RESOLUTION = 0.2;
MAPPING_PLANE_RESOLUTION = 0.4;
HUBER_THRES = 0.1;
%% Open data
dataBag = rosbag('inputData/nsh_indoor_outdoor.bag');
velodyneBag = select(dataBag,'Topic','/velodyne_points');
numMessages = velodyneBag.NumMessages;
%% Run
numMessages = 10;
% preallocation & initialization
systemInited = false;
q_w_curr = [1;0;0;0];
t_w_curr = [0;0;0];
qt_w_curr = [q_w_curr;t_w_curr];
q_last_curr = [1;0;0;0];
t_last_curr = [0;0;0];
qt_last_curr = [q_last_curr;t_last_curr];
pc_cornerLastXyzIdRt = [];
pc_surfLastXyzIdRt = [];
q_w_curr = zeros(numMessages,4);
t_w_curr = zeros(numMessages,3);
trajVisualHolder = gobjects;
mapVisualHolder = gobjects;
pcW = pointCloud(nan(1,3));


timeVar.t_prepare = 0;
timeVar.t_pts = 0;
timeVar.t_opt = 0;
timeVar.t_data = 0;
timeVar.t_solver = 0;

% figure;
% grid on;

timeVar.t_whole = tic;
for epoch = 1:numMessages
    pcRawTemp = readMessages(velodyneBag,epoch);
    pcRawXyz = readXYZ(pcRawTemp{1});
    % scanRegistration
    [pc_sharpXyzIdRt, pc_lessSharpXyzIdRt, pc_flatXyzIdRt, pc_lessFlatXyzIdRt, pcAlignedXyzIdRt, timeVar] ...
        = scanRegistration(pcRawXyz, timeVar, SCAN_LINE, MINIMUM_RANGE);
    % laserOdometry
    [pc_cornerLastXyzIdRt, pc_surfLastXyzIdRt, pcAlignedXyzIdRt, qt_w_curr, qt_last_curr, timeVar, systemInited] ...
        = laserOdometry(pc_sharpXyzIdRt, pc_lessSharpXyzIdRt, pc_flatXyzIdRt, pc_lessFlatXyzIdRt, pc_cornerLastXyzIdRt, pc_surfLastXyzIdRt, pcAlignedXyzIdRt, qt_w_curr, qt_last_curr, timeVar, systemInited, SCAN_LINE, DISTANCE_SQ_THRESHOLD, NEARBY_SCAN, HUBER_THRES);
    q_w_curr(epoch,:) = qt_w_curr(1:4)';
    t_w_curr(epoch,:) = qt_w_curr(5:7)';
    
    %% plot
%     pcW = pcmerge(pcW,pointCloud(rotateframe(quaternion(qt_w_curr(1:4)'),pcAlignedXyzIdRt(:,1:3)) + repmat(qt_w_curr(5:7)',length(pcAlignedXyzIdRt(:,1)),1)),0.2);
%     delete(mapVisualHolder);
%     delete(trajVisualHolder);
%     pcshow(pcW);
%     trajVisualHolder = plot3(t_w_curr(1:epoch,1),t_w_curr(1:epoch,2),t_w_curr(1:epoch,3),'r-','LineWidth',1);
%     hold on;
%     axis equal;
%     xlim([t_w_curr(epoch,1) - 20, t_w_curr(epoch,1) + 20]);
%     ylim([t_w_curr(epoch,2) - 20, t_w_curr(epoch,2) + 20]);
%     title(sprintf('Epoch: %d',epoch));
%     drawnow;
%     pause(0.01);
end
timeVar.t_whole = toc(timeVar.t_whole);

timeVar.t_prepare = timeVar.t_prepare/numMessages;
timeVar.t_pts = timeVar.t_pts/numMessages;
timeVar.t_opt = timeVar.t_opt/(numMessages - 1);
timeVar.t_data = timeVar.t_data/(numMessages - 1);
timeVar.t_solver = timeVar.t_solver/(numMessages - 1);