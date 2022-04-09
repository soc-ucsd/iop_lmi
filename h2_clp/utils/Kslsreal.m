function K = Kslsreal(L,M,R,N,T)
% Compute a state sapce realization of K = L - MR^{-1}N
% Input: 
%       L,M,R,N: FIR factor
%       T      : FIR length
%
% this realization is not minimal
% 
% n: state dimension
% m: input dimension
% p: output dimension

% L \in R^(m,p*(T+1))
% M \in R^(m,n*(T+1))
% R \in R^(n,n*(T+1))
% N \in R^(n,p*(T+1))

[~,p1] = size(L);
[n,~] = size(R);
p      = p1/(T+1);

%% state-space realization for L, M, zR, zN
% L = ss(Zp,hatIp,hatL,L0)
% N = ss(Zp,hatIp,hatN,N0)
hatIp = zeros(T*p,p);  hatIp(1:p,1:p) = eye(p);
Zp   = diag(ones(T-1,1), -1);     % downshift operator
Zp   = kron(Zp,eye(p));
hatL = L(:,p+1:end);
L0   = L(:,1:p);
hatN = N(:,p+1:end);
N0   = N(:,1:p);                   % this one must be zero

% zM = ss(Zn,hatIn,hatM,M1)
% zR = ss(Zn,hatIn,hatR,R1)
hatIn  = zeros((T-1)*n,n); hatIn(1:n,1:n) = eye(n);
Zn = diag(ones(T-2,1), -1);     % downshift operator
Zn = kron(Zn,eye(n));
hatM   = M(:,2*n+1:end);
M1     = M(:,n+1:2*n);
R1     = R(:,n+1:2*n);           % this one must be identity
hatR   = R(:,2*n+1:end);


% state space realization
A = [Zn-hatIn*hatR -hatIn*hatN;
    zeros(p*T,n*(T-1))  Zp];
B = [zeros(n*(T-1),p);hatIp];
C = [hatM-M1*hatR -M1*hatN+hatL];
D = L0;

K = ss(A,B,C,D,[]);

end

