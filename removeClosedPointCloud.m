function ptcloudRawXyz = removeClosedPointCloud(ptcloudRawXyz, MINIMUM_RANGE)
    underThredIdx = sum(ptcloudRawXyz.^2,2) < MINIMUM_RANGE*MINIMUM_RANGE;
    ptcloudRawXyz = ptcloudRawXyz(~underThredIdx,:);
end