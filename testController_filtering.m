% -------------------------------------------------------
%
%    Test Controller design via joint filtering
%
% -------------------------------------------------------



clc

% ----------------------------------
%       Generate ramdon data
% ----------------------------------

% dimension

n = 2;
m = 1;
p = 1;

A = rand(n)*2;
B = rand(n,m);
C = rand(p,n);
D = zeros(p,m);


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

% opts = 1 - joint filtering for performance
% opts = 0 - only find a stabilizing controller

% -------------------------------------------------------
%             Stabilization only 
% -------------------------------------------------------

eps = (0.3)^2;
opts = 0; 
[K,X,Y,mu]     = filterController(Ml,Nl,eps,opts);
[K2,X2,Y2,mu2] = filterController2(Ml,Nl,eps,opts);

% minreal(Y*X^(-1)- K)

G = ss(A,B,C,D,-1);
H = feedback(G,-K);
H2 = feedback(G,-K2);

fprintf('maximal open-loop poles: %6.2f \n',max(abs(pole(G))))
fprintf('Standard Hinf LMI: \n')
fprintf('maximal cloed-loop poles: %6.2f \n',max(abs(pole(H))))
fprintf('Hinf norm - Feasibility: %6.2f \n', hinfnorm([Ml, -Nl]*[X;Y] - eye(size(Ml.C,1))))
fprintf('Hinf norm - Performance: %6.2f \n', hinfnorm([X;Y]*[Ml, Nl] + blkdiag(zeros(size(X.C,1)),eye(size(Y.C,1)))))
fprintf('Extended Hinf LMI: \n')
fprintf('maximal cloed-loop poles: %6.2f \n',max(abs(pole(H2))))
fprintf('Hinf norm - Feasibility: %6.2f \n', hinfnorm([Ml, -Nl]*[X2;Y2] - eye(size(Ml.C,1))))
fprintf('Hinf norm - Performance: %6.2f \n', hinfnorm([X2;Y2]*[Ml, Nl] + blkdiag(zeros(size(X2.C,1)),eye(size(Y2.C,1)))))


fprintf('Paused: type any key to continue \n\n ')

% -------------------------------------------------------
%             optimize performance
% -------------------------------------------------------
opts = 1; 
[K,X,Y,mu]     = filterController(Ml,Nl,eps,opts);
[K2,X2,Y2,mu2] = filterController2(Ml,Nl,eps,opts);

H = feedback(G,-K);
H2 = feedback(G,-K2);

fprintf('Standard Hinf LMI: \n')
fprintf('maximal cloed-loop poles: %6.2f \n',max(abs(pole(H))))
fprintf('Hinf norm - Feasibility: %6.2f \n', hinfnorm([Ml, -Nl]*[X;Y] - eye(size(Ml.C,1))))
fprintf('Hinf norm - Performance: %6.2f %6.2f\n', hinfnorm([X;Y]*[Ml, Nl] + blkdiag(zeros(size(X.C,1)),eye(size(Y.C,1)))),mu)
fprintf('Extended Hinf LMI: \n')
fprintf('maximal cloed-loop poles: %6.2f \n',max(abs(pole(H2))))
fprintf('Hinf norm - Feasibility: %6.2f \n', hinfnorm([Ml, -Nl]*[X2;Y2] - eye(size(Ml.C,1))))
fprintf('Hinf norm - Performance: %6.2f %6.2f \n', hinfnorm([X2;Y2]*[Ml, Nl] + blkdiag(zeros(size(X2.C,1)),eye(size(Y2.C,1)))),mu2)
