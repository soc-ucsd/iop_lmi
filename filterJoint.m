function [Af,Bf,Cf,Df,mu1,mu2] = filterJoint(P,H,mu1,mu2)
% Joint Hinf filtering design
% Input 
%  P =[P1 P2] = [A | B1 B2]
%               [--|------]
%               [C | D1 D2] 
%
%  H = [H1 H2] = [A  | B ]
%               [---|-- ]
%               [C1 | D1] 
%               [C2 | D2]

% Find a common filter F such that 
%      |P1*F - P2|_infty^2 < mu1
%      |F*H1 - H2|_infty^2 < mu2
%
%   or minimize mu1 + mu2 subject to the constriant above
%
% Output: a state-space realization of F
%      F = Cf(zI - Af)^{-1}Bf + Df

if nargin < 3
    flagPerformance = 1;
else
    flagPerformance = 0;
end

% Dimension - right filter
[n,m1] = size(P.B1);       % n - state; m1 - input of P1
[~,m2] = size(P.B2);       % n - state; m2 - input of P2
[p,~]  = size(P.C);        % p - output
 
% Dimension - left filter
[n,m]   = size(H.B);         % n - state; m - input 
[p1,~]  = size(H.C1);        % p1 - output of H1
[p2,~]  = size(H.C2);        % p2 - output of H2
 

% Define variables
X = sdpvar(n);           % symmetric variable - Lyapunov
Z = sdpvar(n);           % symmetric variable - Lyapunov 

Y = sdpvar(n);           % symmetric variable - Lyapunov
%S = sdpvar(n);           % symmetric variable - Lyapunov 

% Filter realization
Q = sdpvar(n,n,'full');
F = sdpvar(n,m2);
L = sdpvar(m1,n);
R = sdpvar(m1,m2);

% Hinf performance
if flagPerformance == 1
    mu1 = sdpvar(1);
    mu2 = sdpvar(1);
end

% Hinf constraint
epsilon = 1e-5;
A = P.A;  B1 = P.B1;  B2 = P.B2;
C = P.C;  D1 = P.D1;  D2 = P.D2;
HinfLMI_right = [X,           Z, A*X+B1*L, A*Z+B1*L, B1*R-B2,      zeros(n,p);
           Z,           Z, Q,         Q,        F,           zeros(n,p);
           (A*X+B1*L)', Q',X,         Z,    zeros(n,m2),     X*C'+L'*D1';
           (A*Z+B1*L)', Q',Z,         Z,    zeros(n,m2),     Z*C'+L'*D1';
           (B1*R-B2)',  F',zeros(m2,n),zeros(m2,n),eye(m2),  R'*D1'-D2';
           zeros(p,n), zeros(p,n),C*X+D1*L,C*Z+D1*L,D1*R-D2, mu1*eye(p)];
   
constraint = [HinfLMI_right - epsilon*eye(4*n+m2+p)>=0];

A = H.A; C1 = H.C1;  C2 = H.C2;
B = H.B; D1 = H.D1;  D2 = H.D2;
HinfLMI_left = [Y,           Z,          Y*A+F*C1,     Q, Y*B+F*D1,    zeros(n,p2);
           Z,           Z,          Z*A+F*C1,     Q, Z*B+F*D1,    zeros(n,p2);
           (Y*A+F*C1)', (Z*A+F*C1)',Y,            Z, zeros(n,m),  C1'*R'-C2';
           Q', Q',                  Z,            Z, zeros(n,m),     L';
           (Y*B+F*D1)', (Z*B+F*D1)',zeros(m,n), zeros(m,n), eye(m), D1'*R'-D2';
           zeros(p2,n), zeros(p2,n),R*C1 - C2,    L,       R*D1-D2, mu2*eye(p2)];
       
constraint = [constraint, HinfLMI_left - epsilon*eye(4*n+m+p2)>=0];

%constraint = [constraint, Z == Z];

% cost - feasibility
if flagPerformance == 0
    cost = 0; %trace(X+Z);
elseif flagPerformance == 1
    cost = mu1 + mu2;
end

% call mosek
ops = sdpsettings('solver','mosek','verbose',0);
optimize(constraint,cost,ops);

% Filter -- state-space realization
vZ = value(Z);
U  = vZ^(0.5);
V  = vZ^(0.5);

Af = U*value(Z)^(-1)*value(Q)*U^(-1);
Bf = U*value(Z)^(-1)*value(F);
Cf = value(L)*U^(-1);
Df = value(R);

% alternative way
% Afa = V^(-1)*value(Q)*vZ^(-1)*V;
% Bfa = V^(-1)*value(F);
% Cfa = value(L)*vZ^(-1)*V;
% Dfa = value(R);

if flagPerformance == 1
    mu1 = value(mu1);
    mu2 = value(mu2);
end

end

