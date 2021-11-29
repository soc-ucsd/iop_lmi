function [F,mu,K,X,Y] = filterleftInverse(P1,P2,P3,Ds)
% Compute a state-space realization
%       min_F |F*P1 - P2|
%              P3*F = I
% Input: state-space model of P1,P2
% Output: state-space model of F


    % state-space model
    A  = P1.A;
    B  = P1.B;
    C1 = P1.C;
    D1 = P1.D;
    C2 = P2.C;
    D2 = P2.D;
    
    B3 = P3.B;
    C3 = P3.C;
    D3 = P3.D;


    % Dimension 
    [n,m]  = size(B3);        % n - state; m - input 
    [p,~]  = size(C3);        % p - output of P

    % Define variables
    P = sdpvar(n);           % symmetric variable - Lyapunov
    K = sdpvar(n);           % symmetric variable - Lyapunov
    
    J = sdpvar(n,n,'full');           % nonsymmetric variable - Lyapunov
    Y = sdpvar(n,n,'full');           % nonsymmetric variable - Lyapunov
    M = sdpvar(n,n,'full');           % nonsymmetric variable - Lyapunov
  
    % Inverse realization
    S = sdpvar(n,n,'full');
    %F = sdpvar(n,p,'full');
    L = sdpvar(m,n,'full');
    R = sdpvar(m,p,'full');

    % cost
    mu = sdpvar(1);
    
    % Constraint 
    epsilon = 1e-5;
    constraint = [C3*S + D3*L == 0, D3*R == eye(p)];

    Hinf = [P,                    J,                 Y*A+B3*R*C1,         A*S+B3*L,      Y*B+B3*R*D1,           zeros(n,size(C2,1));
            J',                   K,                 (Y-M)*A+B3*R*C1,     A*S+B3*L,     (Y-M)*B+B3*R*D1,        zeros(n,size(C2,1));
           (Y*A+B3*R*C1)',      ((Y-M)*A+B3*R*C1)',  Y+Y'-P,              S+(Y-M)'-J,    zeros(n,size(D1,2)),   (R*C1 - C2)';
           (A*S+B3*L)',         (A*S+B3*L)',         S'+(Y-M)-J',         S+S'-K,        zeros(n,size(D1,2)),     L';
           (Y*B+B3*R*D1)',      ((Y-M)*B+B3*R*D1)',  zeros(size(D1,2),n), zeros(size(D1,2),n), eye(size(D1,2)),  (R*D1 - D2)';
           zeros(size(C2,1),n), zeros(size(C2,1),n), R*C1 - C2,            L,             R*D1 - D2,              mu*eye(size(C2,1))];
    
    constraint = [constraint, Hinf-epsilon*eye(4*n+size(C2,1)+size(D1,2)) >=0];

    % cost - feasibility
    cost = mu; 

    % call mosek
    ops = sdpsettings('solver','mosek','verbose',1);
    optimize(constraint,cost,ops);

    % construct a state-space realization
    Af = value(S)^(-1)*(A*value(S)+B3*value(L));
    Bf = value(S)^(-1)*(B3*value(R));
    Cf = value(L);
    Df = value(R);

    F = ss(Af,Bf,Cf,Df,-1);
    mu = value(mu).^(0.5);

    
    K = []; X = []; Y = [];
    
    if nargin >= 4
        % extract controller
        X = ss(Af,Bf,Cf(1:Ds,:),Df(1:Ds,:),-1);
        Y = ss(Af,Bf,Cf(Ds+1:end,:),Df(Ds+1:end,:),-1);

        % controller K = YX^(-1)
        Ak = Af - Bf*X.D^(-1)*X.C;
        Bk = -Bf*X.D^(-1);
        Ck = -Y.C+Y.D*X.D^(-1)*X.C;
        Dk = Y.D*X.D^(-1);
        K = ss(Ak,Bk,Ck,Dk,-1);
    end


end

