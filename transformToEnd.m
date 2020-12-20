function pcXyzIdRTout = transformToEnd(pcXyzIdRTin, q_last_curr, t_last_curr, DISTORTION)
    un_point = transformToStart(pcXyzIdRTin, q_last_curr, t_last_curr, DISTORTION);
    tArr = repmat(t_last_curr',length(pcXyzIdRTin(:,1)),1);
    point_end = rotateframe(conj(quaternion(q_last_curr')),(un_point - tArr));
    pcXyzIdRTout = [point_end,pcXyzIdRTin(:,4:5)];
end