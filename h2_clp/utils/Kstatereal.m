function K = Kstatereal(U,Y,T)
% Compute a state sapce realization of K = UY^(-1)
% Input: 
%       U,Y: FIR factor
%       T  : FIR length
%
% this realization is not minimal
% 
% Check again -- the following realization works now

% U \in R^(m,p*(N+1))
% Y \in R^(p,p*(N+1))

[m,p1] = size(U);
p      = p1/(T+1);

% data in Y, U
hatU = U(:,p+1:end);
U0   = U(:,1:p);
hatY = Y(:,p+1:end);
Y0   = Y(:,1:p);               % this should be identity from IOP

hatI          = zeros(T*p,p);
hatI(1:p,1:p) = eye(p);

Z = diag(ones(T-1,1), -1);     % downshift operator
Z = kron(Z,eye(p));

% state space matrices
A = Z - hatI*hatY;
B = -hatI;
C = U0*hatY - hatU;
D = U0;

K = ss(A,B,C,D,[]);

end

