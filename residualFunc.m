function residues = residualFunc(x, ccp, clpa, clpb, scp, slpj, slpl,slpm, DISTORTION, HUBER_THRES)
    q_last_curr = x(1:4)/norm(x(1:4));
    t_last_curr = x(5:7);
    %% corner points
    q0 = ones(1,'quaternion');
    if(DISTORTION)
        qArr = slerp(q0, quaternion(q_last_curr'),ccp(:,5));
        tArr = ccp(:,5)*t_last_curr';
        clp = rotateframe(qArr,ccp(:,1:3)) + tArr;
    else
        tArr = repmat(t_last_curr',length(ccp(:,1)),1);
        clp = rotateframe(quaternion(q_last_curr'),ccp(:,1:3)) + tArr;
    end
    nu = cross(clp-clpa, clp-clpb, 2);
    de_norm = sqrt(sum((clpa - clpb).^2,2));
    cornerResidues = [nu(:,1)./de_norm, nu(:,2)./de_norm, nu(:,3)./de_norm];
%     cornerResidues = sqrt(sum(nu.^2,2))./de_norm;
    %% surface points
    ljm = cross(slpj - slpl,slpj-slpm,2);
    ljm_norm = sqrt(sum(ljm.^2,2));
    ljm = [ljm(:,1)./ljm_norm, ljm(:,2)./ljm_norm, ljm(:,3)./ljm_norm];
    if(DISTORTION)
        qArr = slerp(q0, quaternion(q_last_curr'),scp(:,5));
        tArr = scp(:,5)*t_last_curr';
        slp = rotateframe(qArr,scp(:,1:3)) + tArr;
    else
        tArr = repmat(t_last_curr',length(scp(:,1)),1);
        slp = rotateframe(quaternion(q_last_curr'),scp(:,1:3)) + tArr;
    end
    surfResidues = sum((slp - slpj).*ljm,2);
    %%
    residues = [reshape(cornerResidues',length(cornerResidues(:,1))*3,1);surfResidues];
%     residues = [cornerResidues;surfResidues];
    %%
    innerSelector = abs(residues)<=HUBER_THRES;
    residues(innerSelector) = sqrt(0.5)*residues(innerSelector);
    residues(~innerSelector) = sqrt(HUBER_THRES*abs(residues(~innerSelector)) - sqrt(0.5)*HUBER_THRES^2);
%     residues(~innerSelector) = 2*sqrt(HUBER_THRES) - 1;
end