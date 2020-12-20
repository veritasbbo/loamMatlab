function pcXyzIdRTout = transformToStart(pcXyzIdRTin, q_last_curr, t_last_curr, DISTORTION)
    q0 = ones(1,'quaternion');
    if(DISTORTION)
        qArr = slerp(q0, quaternion(q_last_curr'),pcXyzIdRTin(:,5));
        tArr = pcXyzIdRTin(:,5)*t_last_curr';
        un_point = rotateframe(qArr,pcXyzIdRTin(:,1:3)) + tArr;
    else
        tArr = repmat(t_last_curr',length(pcXyzIdRTin(:,1)),1);
        un_point = rotateframe(quaternion(q_last_curr'),pcXyzIdRTin(:,1:3)) + tArr;
    end
    pcXyzIdRTout = [un_point,pcXyzIdRTin(:,4:5)];
end