function K = leftInverse(P)
% Compute a state-space realization
%             P*K = I
% Input: state-space model of P
% Output: state-space model of K


    % state-space model
    A = P.A;
    B = P.B;
    C = P.C;
    D = P.D;


    % Dimension 
    [n,m]  = size(B);        % n - state; m - input 
    [p,~]  = size(C);        % p - output of P

    % Define variables
    X = sdpvar(n);           % symmetric variable - Lyapunov

    % Inverse realization
    S = sdpvar(n,n,'full');
    F = sdpvar(n,p,'full');
    L = sdpvar(m,n,'full');
    R = sdpvar(m,p,'full');

    % Constraint 
    epsilon = 1e-5;
    constraint = [B*R == F, C*S + D*L == 0, D*R == eye(p)];

    constraint = [constraint, [X, A*S+B*L; S'*A'+L'*B', S+S'-X]-epsilon*eye(2*n) >=0];

    % cost - feasibility
    cost = 0; %trace(X+Z);

    % call mosek
    ops = sdpsettings('solver','mosek','verbose',1);
    optimize(constraint,cost,ops);

    % construct a state-space realization
    Ak = value(S)^(-1)*(A*value(S)+B*value(L));
    Bk = value(S)^(-1)*(value(F));
    Ck = value(L);
    Dk = value(R);

    K = ss(Ak,Bk,Ck,Dk,-1);

end

