function [K,X,Y] = filterControllerDec(Ml,Nl,eps,opts)
% Decentralized controller design via right filtering
% Input: 
%       Left Co-prime factor - G = Ml^(-1)Nl
%                     [Ml Nl] = [A | Bm Bn]
%                               [---------]
%                               [C | Dm Dn]
%       eps <= 1       : enforce stability 
%       opts.decen    0: stablization only 
%                     1: Decentralized stabilization
% Output

% Used the standard Hinf LMI

if nargin < 4
    opts.decen = 0;      % stabilization or performance
end

n = size(Ml.A,1);

% Right filtering - enforce stability
P.A  = Ml.A; P.B1  = [Ml.B, -Nl.B];
P.C  = Ml.C; P.B2 = zeros(n,size(Ml.C,1));
P.D1 = [Ml.D, -Nl.D];
P.D2 = eye(size(Ml.C,1));


% Dimension - right filter
[n,m1] = size(P.B1);       % n - state; m1 - input of P1
[~,m2] = size(P.B2);       % n - state; m2 - input of P2
[p,~]  = size(P.C);        % p - output
 

if opts.decen == 0 % not decentralized
    % Define variables
    X = sdpvar(n);           % symmetric variable - Lyapunov
    Z = sdpvar(n);           % symmetric variable - Lyapunov 

    % Filter realization
    Q = sdpvar(n,n,'full');
    F = sdpvar(n,m2,'full');
    L = sdpvar(m1,n,'full');
    R = sdpvar(m1,m2,'full');

else
    %X = [];           % symmetric variable - Lyapunov
    X = sdpvar(n);
    Z = [];           % symmetric variable - Lyapunov 

    % Filter realization
    Q = [];
    F = [];
    L1 = []; L2 = [];
    R1 = []; R2 = [];
    
    for i = 1: opts.agent
        %X = blkdiag(X,sdpvar(opts.n));
        Z = blkdiag(Z,sdpvar(opts.n));

        % Filter realization
        Q  = blkdiag(Q,sdpvar(opts.n,opts.n,'full'));
        F  = blkdiag(F,sdpvar(opts.n,opts.p,'full'));
        L1 = blkdiag(L1,sdpvar(opts.p,opts.n,'full'));
        R1 = blkdiag(R1,sdpvar(opts.p,opts.p,'full'));
        L2 = blkdiag(L2,sdpvar(opts.m,opts.n,'full'));
        R2 = blkdiag(R2,sdpvar(opts.m,opts.p,'full'));
    end
    L = [L1;L2];
    R = [R1;R2];
    
end


% enforce stability -- eps < 1
epsilon = 1e-5;
A = P.A;  B1 = P.B1;  B2 = P.B2;
C = P.C;  D1 = P.D1;  D2 = P.D2;
HinfLMI_right = [X,           Z, A*X+B1*L, A*Z+B1*L, B1*R-B2,      zeros(n,p);
           Z,           Z, Q,         Q,        F,           zeros(n,p);
           (A*X+B1*L)', Q',X,         Z,    zeros(n,m2),     X*C'+L'*D1';
           (A*Z+B1*L)', Q',Z,         Z,    zeros(n,m2),     Z*C'+L'*D1';
           (B1*R-B2)',  F',zeros(m2,n),zeros(m2,n),eye(m2),  R'*D1'-D2';
           zeros(p,n), zeros(p,n),C*X+D1*L,C*Z+D1*L,D1*R-D2, eps*eye(p)];
   
%constraint = [HinfLMI_right - epsilon*eye(4*n+m2+p)>=0];
constraint = [HinfLMI_right >=0, X - Z >= epsilon*eye(n)];

if opts.cost == 1       % stabilization only 
    cost = 0;
elseif opts.cost == 2   % regularize state-space realization
    cost = norm(Q,inf) + norm(F,inf) + norm(L,inf) + norm(R,inf);
end
    
bigM = 1e8;
constraint = [constraint, norm([Q, F;L R],2) <= bigM];
%                           norm(F,1) <= bigM,
%                           norm(L,1) <= bigM,
%                           norm(R,1) <= bigM];
                    
constraint = [constraint, norm(X,2) <= bigM];
% call mosek
ops = sdpsettings('solver','mosek','verbose',1);
optimize(constraint,cost,ops);

% Filter -- state-space realization
vZ = value(Z);
U  = vZ^(0.5);
V  = vZ^(0.5);

Af = U*value(Z)^(-1)*value(Q)*U^(-1);
Bf = U*value(Z)^(-1)*value(F);
Cf = value(L)*U^(-1);
Df = value(R);

X = ss(Af,Bf,Cf(1:size(Ml.D,1),:),Df(1:size(Ml.D,1),:),-1);
Y = ss(Af,Bf,Cf(size(Ml.D,1)+1:end,:),Df(size(Ml.D,1)+1:end,:),-1);

% controller realization K = YX^(-1)
Ak = Af - Bf*X.D^(-1)*X.C;
Bk = -Bf*X.D^(-1);
Ck = -Y.C+Y.D*X.D^(-1)*X.C;
Dk = Y.D*X.D^(-1);
K  = ss(Ak,Bk,Ck,Dk,-1);


end



