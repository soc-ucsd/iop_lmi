%
% Controller design via joint filtering

clc

% dimension
n = 2;
m = 1;
p = 1;

A = rand(n)*2;
B = rand(n,m);
C = rand(p,n);
D = rand(p,m);


% Step 1: Left co-prime factorization
poles = rand(n,1) - 0.5;
K = place(A,B,poles);
F = place(A',C',poles);
L = -F';

Ml = ss(A + L*C, L, C, eye(p),-1);  % G = Ml^(-1)*Nl
Nl = ss(A + L*C, B+L*D, C, D,-1);

Mr = ss(A - B*K, B, -K, eye(m),-1); % G = Nr*Mr^(-1)
Nr = ss(A - B*K, B,C-D*K, D,-1);



% Step 2: enforce stability
P = ss(Ml.A,[Ml.B, -Nl.B],Ml.C,[Ml.D, -Nl.D],-1);
K0 = leftInverse(P);

minreal(P*K0-eye(p))

% G = ss(A,B,C,D,-1);
% H = feedback(G,-K0);
% pole(G)
% pole(H)
% 
% eps = 0.2;
% [K,X,Y,mu] = filterController(Ml,Nl,eps,0);
% 
% minreal(Y*X^(-1)- K)


