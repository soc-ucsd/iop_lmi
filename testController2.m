% -------------------------------------------------------------------------
%    Test Controller design via filtering with equality constraints
% -------------------------------------------------------------------------

clc

% dimension
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


P1 = ss(Ml.A,[Ml.B, Nl.B],Ml.C,[Ml.D, Nl.D],-1); 
P2 = ss([],[],[],[],-1);
P2.A = Ml.A; P2.B  = [Ml.B, Nl.B];
P2.C = zeros(size(Ml.B,2)+size(Nl.B,2),n);
P2.D = zeros(size(Ml.B,2)+size(Nl.B,2));
P2.D(size(Ml.B,2)+1:end,size(Ml.B,2)+1:end) = -eye(size(Nl.B,2));

P3 = ss(Ml.A,[Ml.B, -Nl.B],Ml.C,[Ml.D, -Nl.D],-1);

%  minreal(((Nr*Mr^(-1))) - (ss(A,B,C,D)),1e-6)
%  minreal(((Ml^(-1)*Nl)) - (ss(A,B,C,D)),1e-6)

[F,mu,K,X,Y] = filterleftInverse(P1,P2,P3,size(Ml.D,1));

eps = (0.2)^2;
[K1,X1,Y1,mu1] = filterController(Ml,Nl,eps,1);
[K2,X2,Y2,mu2] = filterController2(Ml,Nl,eps,1);


fprintf('Filter with inverse, Joint filtering (LMI), Joint fitering (extended LMI)\n\n')
fprintf('Feasibility - Hinf norm: %6.4f, %6.4f, %6.4f\n',hinfnorm(P3*F - eye(p)),hinfnorm(P3*[X1;Y1] - eye(p)),hinfnorm(P3*[X2;Y2] - eye(p)))
fprintf('Performance - mu: %6.4f, %6.4f, %6.4f\n',mu,mu1,mu2)


% %% Hinf control toolbox
% % x[t+1] = Ax[t] + Bu[t] + B*du
% % z[t]   = [y[t];
% %           u[t]]
% % y[t]   = Cx[t] + Du[t] + dy
% 
% % w[t] = [du;dy]
% 
% % Generalized state-space model
% % x[t+1] = A x[t] + B1 w[t] + B2  u[t]
% % z[t]   = C1x[t] + D11w[t] + D12 u[t]
% % y[t]   = C2x[t] + D21w[t] + D22 u[t]
% 
% B1  = [B zeros(n,p)];      B2 = B;
% D21 = [zeros(p,m) eye(p)]; D22 = D;
% C2  = C;
% C1  = [C2;
%        zeros(m,n)];
% D11 = [D21;
%        zeros(m,p+m)];
% D12 = [zeros(p,m);
%         eye(m)];
% %  generalized state-space model
% Bg = [B1, B2];
% Cg = [C1;C2];
% Dg = [D11 D12;
%       D21 D22];
% 
% P = ss(A,Bg,Cg,Dg,-1);            % discrete time model -- transfer matrices
% [K3,CL,gamma,info] = hinfsyn(P,p,m);  % this norm is not consistent with the results from CLP. 
% 
% [mu,mu1,mu2,gamma]
% 
% 
% [hinfnorm([eye(p), -G; -K, eye(m)]^(-1));
%  hinfnorm([eye(p), -G; -K1, eye(m)]^(-1));
%  hinfnorm([eye(p), -G; -K2, eye(m)]^(-1));
%  hinfnorm([eye(p), -G; -K3, eye(m)]^(-1))]