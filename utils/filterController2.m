function [K,X,Y,mu] = filterController2(Ml,Nl,eps,opts)
% Hinf controller design via joint filtering
% Input: 
%       Left Co-prime factor - G = Ml^(-1)Nl
%                     [Ml Nl] = [A | Bm Bn]
%                               [---------]
%                               [C | Dm Dn]
%       eps < 1 - enforce stability 
%  
% Output

% used the extended Hinf inequality in Theorem 2 from
%        De Oliveira, M. C., Geromel, J. C., & Bernussou, J. (2002). 
%          Extended H 2 and H norm characterizations and controller parametrizations 
%          for discrete-time systems. International journal of control, 75(9), 666-679.


if nargin < 4
    opts = 1; % performance
end

n = size(Ml.A,1);

% Right filtering - enforce stability
Pr.A  = Ml.A; Pr.B1  = [Ml.B, -Nl.B];
Pr.C  = Ml.C; Pr.B2 = zeros(n,size(Ml.C,1));
Pr.D1 = [Ml.D, -Nl.D];
Pr.D2 = eye(size(Ml.C,1));


% Left filtering - encode Hinf performance
H.A  = Ml.A; H.B  = [Ml.B, Nl.B];
H.C1 = Ml.C; H.D1 = [Ml.D, Nl.D];
H.C2 = zeros(size(Ml.B,2)+size(Nl.B,2),n);
H.D2 = zeros(size(Ml.B,2)+size(Nl.B,2));
H.D2(size(Ml.B,2)+1:end,size(Ml.B,2)+1:end) = -eye(size(Nl.B,2));


% Dimension - right filter
[n,m1] = size(Pr.B1);       % n - state; m1 - input of P1
[~,m2] = size(Pr.B2);       % n - state; m2 - input of P2
[p,~]  = size(Pr.C);        % p - output
 
% Dimension - left filter
[n,m]   = size(H.B);         % n - state; m - input 
[p1,~]  = size(H.C1);        % p1 - output of H1
[p2,~]  = size(H.C2);        % p2 - output of H2
 

if opts == 1      % Hinf performance design
    
elseif opts == 0  % enforce stability only
    
end

% Define variables
E = sdpvar(n);            % symmetric variable - Lyapunov
Hr = sdpvar(n);           % symmetric variable - Lyapunov 

X = sdpvar(n,n,'full');           % nonsymmetric variable
Z = sdpvar(n,n,'full');           % nonsymmetric variable  
N = sdpvar(n,n,'full');           % nonsymmetric variable
G = sdpvar(n,n,'full');           % nonsymmetric variable  
    
% Left Hinf
P = sdpvar(n);
K = sdpvar(n);

Y = sdpvar(n,n,'full');           % nonsymmetric variable
J = sdpvar(n,n,'full');           % nonsymmetric variable
M = sdpvar(n,n,'full');           % nonsymmetric variable  


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
A = Pr.A;  B1 = Pr.B1;  B2 = Pr.B2;
C = Pr.C;  D1 = Pr.D1;  D2 = Pr.D2;
    HinfLMI_right = [E,              G, A*X+B1*L, A*(X-N)+B1*L, B1*R-B2,      zeros(n,p);
               G',             Hr, Q,         Q,        F,              zeros(n,p);
              (A*X+B1*L)',     Q',X+X'-E, X-N+Z'-G,    zeros(n,m2),     X'*C'+L'*D1';
              (A*(X-N)+B1*L)', Q',(X-N+Z'-G)', Z+Z'-Hr,  zeros(n,m2),  (X-N)'*C'+L'*D1';
              (B1*R-B2)',      F',zeros(m2,n),zeros(m2,n),eye(m2),     R'*D1'-D2';
              zeros(p,n), zeros(p,n), C*X+D1*L, C*(X-N)+D1*L,   D1*R-D2,  eps*eye(p)];

    %constraint = [HinfLMI_right - epsilon*eye(4*n+m2+p)>=0];
    constraint = [HinfLMI_right >=0, [E, G; G', Hr]-epsilon*eye(2*n)>=0];

if opts == 1 % Hinf performance
    % Hinf performance optimization
    A = H.A; C1 = H.C1;  C2 = H.C2;
    B = H.B; D1 = H.D1;  D2 = H.D2;
    HinfLMI_left = [P,           J,          Y*A+F*C1,     Q, Y*B+F*D1,         zeros(n,p2);
               J',           K,     (Y-M)*A+F*C1,     Q, (Y-M)*B+F*D1,     zeros(n,p2);
              (Y*A+F*C1)', ((Y-M)*A+F*C1)',Y+Y'-P, Z+(Y-M)'-J, zeros(n,m), C1'*R'-C2';
               Q', Q',         Z'+(Y-M)-J',         Z'+Z-K, zeros(n,m),     L';
               (Y*B+F*D1)', ((Y-M)*B+F*D1)',zeros(m,n), zeros(m,n), eye(m), D1'*R'-D2';
               zeros(p2,n), zeros(p2,n),R*C1 - C2,    L,       R*D1-D2, mu*eye(p2)];
    %constraint = [constraint, HinfLMI_left - epsilon*eye(4*n+m+p2)>=0];   
    
    
    constraint = [constraint, HinfLMI_left >=0, [P, J; J', K]-epsilon*eye(2*n)>=0];
    % cost - feasibility
    cost = mu;
end

%bigM = 1e6;
bigM = 1e8;
constraint = [constraint, norm([Q, F;L R],2) <= bigM];
% constraint = [constraint, norm(Q,1) <= bigM,
%                           norm(F,1) <= bigM,
%                           norm(L,1) <= bigM,
%                           norm(R,1) <= bigM];
  

constraint = [constraint, norm([X, X-N;Z Z],2) <= bigM];
% constraint = [constraint, norm(E,1) <= bigM,
%                           norm(Hr,1) <= bigM,
%                           norm(X,1) <= bigM,
%                           norm(Z,1) <= bigM,
%                           norm(N,1) <= bigM,
%                           norm(G,1) <= bigM,
%                           norm(E,1) <= bigM,
%                           norm(P,1) <= bigM,
%                           norm(K,1) <= bigM,
%                           norm(Y,1) <= bigM,
%                           norm(J,1) <= bigM,
%                           norm(M,1) <= bigM];

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



