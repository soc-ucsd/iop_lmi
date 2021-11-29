function [K,X,Y,mu] = filterController(Ml,Nl,eps,opts)
% Hinf controller design via joint filtering
% Input: 
%       Left Co-prime factor - G = Ml^(-1)Nl
%                     [Ml Nl] = [A | Bm Bn]
%                               [---------]
%                               [C | Dm Dn]
%       eps < 1 - enforce stability 
%  
% Output

if nargin < 4
    opts = 1; % performance
end

n = size(Ml.A,1);

% Right filtering - enforce stability
P.A  = Ml.A; P.B1  = [Ml.B, -Nl.B];
P.C  = Ml.C; P.B2 = zeros(n,size(Ml.C,1));
P.D1 = [Ml.D, -Nl.D];
P.D2 = eye(size(Ml.C,1));


% Left filtering - encode Hinf performance
H.A  = Ml.A; H.B  = [Ml.B, Nl.B];
H.C1 = Ml.C; H.D1 = [Ml.D, Nl.D];
H.C2 = zeros(size(Ml.B,2)+size(Nl.B,2),n);
H.D2 = zeros(size(Ml.B,2)+size(Nl.B,2));
H.D2(size(Ml.B,2)+1:end,size(Ml.B,2)+1:end) = -eye(size(Nl.B,2));


% Dimension - right filter
[n,m1] = size(P.B1);       % n - state; m1 - input of P1
[~,m2] = size(P.B2);       % n - state; m2 - input of P2
[p,~]  = size(P.C);        % p - output
 
% Dimension - left filter
[n,m]   = size(H.B);         % n - state; m - input 
[p1,~]  = size(H.C1);        % p1 - output of H1
[p2,~]  = size(H.C2);        % p2 - output of H2
 

if opts == 1      % Hinf performance design
    
elseif opts == 0  % enforce stability only
    
end

% Define variables
X = sdpvar(n);           % symmetric variable - Lyapunov
Z = sdpvar(n);           % symmetric variable - Lyapunov 

Y = sdpvar(n);           % symmetric variable - Lyapunov
%S = sdpvar(n);          % symmetric variable - Lyapunov 

% Filter realization
Q = sdpvar(n,n,'full');
F = sdpvar(n,m2,'full');
L = sdpvar(m1,n,'full');
R = sdpvar(m1,m2,'full');

% Hinf performance
mu = sdpvar(1);
cost = 0;

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
   
constraint = [HinfLMI_right - epsilon*eye(4*n+m2+p)>=0];

if opts == 1 % Hinf performance
    % Hinf performance optimization
    A = H.A; C1 = H.C1;  C2 = H.C2;
    B = H.B; D1 = H.D1;  D2 = H.D2;
    HinfLMI_left = [Y,           Z,          Y*A+F*C1,     Q, Y*B+F*D1,    zeros(n,p2);
               Z,           Z,          Z*A+F*C1,     Q, Z*B+F*D1,    zeros(n,p2);
               (Y*A+F*C1)', (Z*A+F*C1)',Y,            Z, zeros(n,m),  C1'*R'-C2';
               Q', Q',                  Z,            Z, zeros(n,m),     L';
               (Y*B+F*D1)', (Z*B+F*D1)',zeros(m,n), zeros(m,n), eye(m), D1'*R'-D2';
               zeros(p2,n), zeros(p2,n),R*C1 - C2,    L,       R*D1-D2, mu*eye(p2)];

    constraint = [constraint, HinfLMI_left - epsilon*eye(4*n+m+p2)>=0];
    
    % cost - feasibility
    cost = mu;
end

M = 1e4;
constraint = [constraint, norm(Q,1) <= M,
                          norm(F,1) <= M,
                          norm(L,1) <= M,
                          norm(R,1) <= M];

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
mu = value(mu).^(0.5);

% controller K = YX^(-1)
Ak = Af - Bf*X.D^(-1)*X.C;
Bk = -Bf*X.D^(-1);
Ck = -Y.C+Y.D*X.D^(-1)*X.C;
Dk = Y.D*X.D^(-1);
K = ss(Ak,Bk,Ck,Dk,-1);


% alternative way
% Afa = V^(-1)*value(Q)*vZ^(-1)*V;
% Bfa = V^(-1)*value(F);
% Cfa = value(L)*vZ^(-1)*V;
% Dfa = value(R);


end



