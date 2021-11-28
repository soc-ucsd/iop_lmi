function [Af,Bf,Cf,Df] = filterLeft(A,B,C1,C2,D1,D2,mu)
% Right Hinf filtering design
% Input 
%     [H1 H2] = [A  | B ]
%               [---|-- ]
%               [C1 | D1] 
%               [C2 | D2]
% Find a filter F such that 
%      |F*H1 - H2|_infty^2 < mu
%
% Output: a state-space realization of F
%      F = Cf(zI - Af)^{-1}Bf + Df


% Dimension 
[n,m]   = size(B);         % n - state; m - input 
[p1,~]  = size(C1);        % p1 - output of P1
[p2,~]  = size(C2);        % p2 - output of P2
 
% Define variables
Y = sdpvar(n);           % symmetric variable - Lyapunov
S = sdpvar(n);           % symmetric variable - Lyapunov 

% Filter realization
Q = sdpvar(n,n,'full');
F = sdpvar(n,p1);
L = sdpvar(p2,n);
R = sdpvar(p2,p1);

% Hinf constraint
epsilon = 1e-5;
HinfLMI = [Y,           S,          Y*A+F*C1,     Q, Y*B+F*D1,    zeros(n,p2);
           S,           S,          S*A+F*C1,     Q, S*B+F*D1,    zeros(n,p2);
           (Y*A+F*C1)', (S*A+F*C1)',Y,            S, zeros(n,m),  C1'*R'-C2';
           Q', Q',                  S,            S, zeros(n,m),     L';
           (Y*B+F*D1)', (S*B+F*D1)',zeros(m,n), zeros(m,n), eye(m), D1'*R'-D2';
           zeros(p2,n), zeros(p2,n),R*C1 - C2,    L,       R*D1-D2, mu*eye(p2)];
       
constraint = [HinfLMI - epsilon*eye(4*n+m+p2)>=0];

% cost - feasibility
cost = 0; %trace(X+Z);

% call mosek
ops = sdpsettings('solver','mosek','verbose',0);
optimize(constraint,cost,ops);

% Filter -- state-space realization
Af = value(Q)*value(S)^(-1);
Bf = value(F);
Cf = value(L)*value(S)^(-1);
Df = value(R);


end

