function [Af,Bf,Cf,Df] = filterRight(A,B1,B2,C,D1,D2,mu,opts)
% Right Hinf filtering design
% Input 
%     [P1 P2] = [A | B1 B2]
%               [--|------]
%               [C | D1 D2] 
% Find a filter F such that 
%      |P1F - P2|_infty^2 < mu
%
% Output: a state-space realization of F
%      F = Cf(zI - Af)^{-1}Bf + Df

% opts: 
%    1 - standard Hinf inequality 
%    2 - new Hinf inequality in Theorem 2
%        De Oliveira, M. C., Geromel, J. C., & Bernussou, J. (2002). 
%          Extended H 2 and H norm characterizations and controller parametrizations 
%          for discrete-time systems. International journal of control, 75(9), 666-679.


if nargin < 8
    opts = 1;
end

% Dimension 
[n,m1] = size(B1);       % n - state; m1 - input of P1
[~,m2] = size(B2);       % n - state; m2 - input of P2
[p,~]  = size(C);        % p - output
 

% Filter realization
Q = sdpvar(n,n,'full');
F = sdpvar(n,m2,'full');
L = sdpvar(m1,n,'full');
R = sdpvar(m1,m2,'full');

% Hinf constraint

if opts == 1
    % Define variables
    X = sdpvar(n);           % symmetric variable - Lyapunov
    Z = sdpvar(n);           % symmetric variable - Lyapunov 

    epsilon = 1e-5;
    HinfLMI = [X,           Z, A*X+B1*L, A*Z+B1*L, B1*R-B2,      zeros(n,p);
               Z,           Z, Q,         Q,        F,           zeros(n,p);
               (A*X+B1*L)', Q',X,         Z,    zeros(n,m2),     X*C'+L'*D1';
               (A*Z+B1*L)', Q',Z,         Z,    zeros(n,m2),     Z*C'+L'*D1';
               (B1*R-B2)',  F',zeros(m2,n),zeros(m2,n),eye(m2),  R'*D1'-D2';
               zeros(p,n), zeros(p,n),C*X+D1*L,C*Z+D1*L,D1*R-D2, mu*eye(p)];
elseif opts == 2
    % Define variables
    E = sdpvar(n);           % symmetric variable - Lyapunov
    H = sdpvar(n);           % symmetric variable - Lyapunov 
    
    X = sdpvar(n,n,'full');           % nonsymmetric variable
    Z = sdpvar(n,n,'full');           % nonsymmetric variable  
    N = sdpvar(n,n,'full');           % nonsymmetric variable
    G = sdpvar(n,n,'full');           % nonsymmetric variable  


    epsilon = 1e-5;
    HinfLMI = [E,              G, A*X+B1*L, A*(X-N)+B1*L, B1*R-B2,      zeros(n,p);
               G',             H, Q,         Q,        F,              zeros(n,p);
              (A*X+B1*L)',     Q',X+X'-E, X-N+Z'-G,    zeros(n,m2),     X'*C'+L'*D1';
              (A*(X-N)+B1*L)', Q',(X-N+Z'-G)', Z+Z'-H,  zeros(n,m2),  (X-N)'*C'+L'*D1';
              (B1*R-B2)',      F',zeros(m2,n),zeros(m2,n),eye(m2),     R'*D1'-D2';
              zeros(p,n), zeros(p,n), C*X+D1*L, C*(X-N)+D1*L,   D1*R-D2,  mu*eye(p)];
end
       
constraint = [HinfLMI - epsilon*eye(4*n+m2+p)>=0];

% cost - feasibility
cost = 0; %trace(X+Z);

% call mosek
ops = sdpsettings('solver','mosek','verbose',0);
optimize(constraint,cost,ops);

% Filter -- state-space realization
Af = value(Z)^(-1)*value(Q);
Bf = value(Z)^(-1)*value(F);
Cf = value(L);
Df = value(R);


end

