function CL = closedloop(G,K)
% Compute the closed-loop responses
% G: plant, K: controller

A  = G.A;B  = G.B;C  = G.C;D  = G.D;
Ak = K.A;Bk = K.B;Ck = K.C;Dk = K.D;

n  = size(A,1);   % system state
m  = size(B,2);   % system input
p  = size(C,1);   % system output

nk = size(Ak,1);  % controller state

Dinv = [eye(m) -Dk;
        -D     eye(p)]^(-1);
hatC = [zeros(size(Ck,1),n) Ck;
        C      zeros(p,nk)];

Ac = blkdiag(A,Ak) + blkdiag(B,Bk)*Dinv*hatC;
Bc = blkdiag(B,Bk)*Dinv*[zeros(m,p);eye(p)];
Cc = Dinv*hatC;
Dc = Dinv*[zeros(m,p);eye(p)];

CL = ss(Ac,Bc,Cc,Dc,[]);   % closed-loop system



end

