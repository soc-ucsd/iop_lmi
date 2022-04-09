function [info,K] = mixIIpost(XYv,UYv,XUv,UUv,T,m,p,info)
% iop - post processing

        info.var.XY  = (value(XYv));
        info.var.UY  = (value(UYv));
        info.var.XU  = (value(XUv));
        info.var.UU  = (value(UUv));
        
        
        %% state-space realization 
        
        % data in Y, U
        hatU = info.var.UY(:,p+1:end);
        U0   = info.var.UY(:,1:p);
        hatY = info.var.UU(:,m+1:end);
        Y0   = info.var.UU(:,1:m);               % this should be identity from IOP

        hatIp          = zeros(T*p,p);
        hatIp(1:p,1:p) = eye(p);

        Zp = diag(ones(T-1,1), -1);     % downshift operator
        Zp = kron(Zp,eye(p));

        hatIm          = zeros(T*m,m);
        hatIm(1:m,1:m) = eye(m);

        Zm = diag(ones(T-1,1), -1);     % downshift operator
        Zm = kron(Zm,eye(m));
        
        % state space matrices
        A = [Zm-hatIm*hatY -hatIm*hatU;
             zeros(p*T,m*T) Zp];
        B = [-hatIm*U0; hatIp];
        C = [hatY hatU];
        D = U0;

        K = ss(A,B,C,D,[]);
        info.Ks     = K;
        
 
end

