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


%  minreal(((Nr*Mr^(-1))) - (ss(A,B,C,D)),1e-6)
%  minreal(((Ml^(-1)*Nl)) - (ss(A,B,C,D)),1e-6)

eps = (0.1)^2;
[K,X,Y,mu] = filterController(Ml,Nl,eps,1);
[K2,X2,Y2,mu2] = filterController2(Ml,Nl,eps,1);

%minreal(Y*X^(-1)- K)

G = ss(A,B,C,D,-1);
H = feedback(G,-K);
H2 = feedback(G,-K2);
[max(abs(pole(G))),max(abs(pole(H))),max(abs(pole(H2)))]

[hinfnorm([Ml, -Nl]*[X;Y] - eye(size(Ml.C,1))), hinfnorm([X;Y]*[Ml, Nl] + blkdiag(zeros(size(X.C,1)),eye(size(Y.C,1)))),mu]

[hinfnorm([Ml, -Nl]*[X2;Y2] - eye(size(Ml.C,1))), hinfnorm([X2;Y2]*[Ml, Nl] + blkdiag(zeros(size(X2.C,1)),eye(size(Y2.C,1)))),mu2]
% [mu,(1+sqrt(eps)/(1-sqrt(eps)))*mu,hinfnorm([eye(p) -G; -K eye(m)]^(-1))]
