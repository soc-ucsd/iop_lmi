function [Af,Bf,Cf,Df,mu1,mu2] = filterJoint(P,H,opts,mu1,mu2)
% Joint Hinf filtering design
% Input 
%  P = [P1 P2] = [A | B1 B2]
%                [--|------]
%                [C | D1 D2] 
%
%  H = [H1 H2] = [A  | B ]
%                [---|-- ]
%                [C1 | D1] 
%                [C2 | D2]

% Find a common filter F such that 
%      |P1*F - P2|_infty^2 < mu1
%      |F*H1 - H2|_infty^2 < mu2
%
%   or minimize mu1 + mu2 subject to the constriants above
%
% Output: a state-space realization of F
%      F = Cf(zI - Af)^{-1}Bf + Df

% opts: 
%    1 - standard Hinf inequality 
%    2 - new Hinf inequality in Theorem 2
%        De Oliveira, M. C., Geromel, J. C., & Bernussou, J. (2002). 
%          Extended H 2 and H norm characterizations and controller parametrizations 
%          for discrete-time systems. International journal of control, 75(9), 666-679.

if nargin < 4
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
 

% Filter realization
Q = sdpvar(n,n,'full');
F = sdpvar(n,m2,'full');
L = sdpvar(m1,n,'full');
R = sdpvar(m1,m2,'full');

% Hinf performance 
if flagPerformance == 1
    mu1 = sdpvar(1);
    mu2 = sdpvar(1);
end


if opts == 1  % Standard Hinf LMI
    X = sdpvar(n);           % symmetric variable - Lyapunov
    Z = sdpvar(n);           % symmetric variable - Lyapunov
    Y = sdpvar(n);           % symmetric variable - Lyapunov

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

elseif opts == 2 % new Hinf LMI
   
    % Right Hinf 

    A = P.A;  B1 = P.B1;  B2 = P.B2;
    C = P.C;  D1 = P.D1;  D2 = P.D2;
    E = sdpvar(n);            % symmetric variable - Lyapunov
    Hr = sdpvar(n);           % symmetric variable - Lyapunov 
    
    X = sdpvar(n,n,'full');           % nonsymmetric variable
    Z = sdpvar(n,n,'full');           % nonsymmetric variable  
    N = sdpvar(n,n,'full');           % nonsymmetric variable
    G = sdpvar(n,n,'full');           % nonsymmetric variable  

    epsilon = 1e-5;
    HinfLMI_right = [E,              G, A*X+B1*L, A*(X-N)+B1*L, B1*R-B2,      zeros(n,p);
               G',             Hr, Q,         Q,        F,              zeros(n,p);
              (A*X+B1*L)',     Q',X+X'-E, X-N+Z'-G,    zeros(n,m2),     X'*C'+L'*D1';
              (A*(X-N)+B1*L)', Q',(X-N+Z'-G)', Z+Z'-Hr,  zeros(n,m2),  (X-N)'*C'+L'*D1';
              (B1*R-B2)',      F',zeros(m2,n),zeros(m2,n),eye(m2),     R'*D1'-D2';
              zeros(p,n), zeros(p,n), C*X+D1*L, C*(X-N)+D1*L,   D1*R-D2,  mu1*eye(p)];

    constraint = [HinfLMI_right - epsilon*eye(4*n+m2+p)>=0];

  
    % Left Hinf
    P = sdpvar(n);
    K = sdpvar(n);
    
    Y = sdpvar(n,n,'full');           % nonsymmetric variable
    J = sdpvar(n,n,'full');           % nonsymmetric variable
    M = sdpvar(n,n,'full');           % nonsymmetric variable  


    % Hinf constraint
    A = H.A; C1 = H.C1;  C2 = H.C2;
    B = H.B; D1 = H.D1;  D2 = H.D2;
    HinfLMI_left = [P,           J,          Y*A+F*C1,     Q, Y*B+F*D1,         zeros(n,p2);
               J',           K,     (Y-M)*A+F*C1,     Q, (Y-M)*B+F*D1,     zeros(n,p2);
              (Y*A+F*C1)', ((Y-M)*A+F*C1)',Y+Y'-P, Z+(Y-M)'-J, zeros(n,m), C1'*R'-C2';
               Q', Q',         Z'+(Y-M)-J',         Z'+Z-K, zeros(n,m),     L';
               (Y*B+F*D1)', ((Y-M)*B+F*D1)',zeros(m,n), zeros(m,n), eye(m), D1'*R'-D2';
               zeros(p2,n), zeros(p2,n),R*C1 - C2,    L,       R*D1-D2, mu2*eye(p2)];
    constraint = [constraint, HinfLMI_left - epsilon*eye(4*n+m+p2)>=0];   
end
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

% alternative way -- equivalent to the formula above
% Afa = V^(-1)*value(Q)*vZ^(-1)*V;
% Bfa = V^(-1)*value(F);
% Cfa = value(L)*vZ^(-1)*V;
% Dfa = value(R);

if flagPerformance == 1
    mu1 = value(mu1);
    mu2 = value(mu2);
end

end

