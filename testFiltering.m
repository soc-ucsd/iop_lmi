
% ----------------------------------------------
%    Test LMI formulations for filter design
% ----------------------------------------------

clc;clear

%% Right filtering
n  = 6; 
m1 = 2;
m2 = 1; 
p  = 3;

Ar = rand(n);
Ar = Ar/max(abs(eig(Ar)))/1.1;  % stable A

B1 = rand(n,m1);
B2 = rand(n,m2);
C  = rand(p,n);
D1r = zeros(p,m1);
D2r = zeros(p,m2);

% construct the problem data P
% P = [P1 P2] = [A | B1 B2]
%               [--|------]
%               [C | D1 D2] 
% |P1F - P2|_infty^2 < mu

P.A = Ar;P.B1 = B1;  P.B2 = B2;
P.C = C; P.D1 = D1r; P.D2 = D2r;

fprintf('Right filtering with standard and extended Hinf LMIs \n ')

mu1  = 1;  % feasibility test
opts = 1;  % 1 - standard Hinf LMI; 2 - extended Hinf inequality
[Af,Bf,Cf,Df] = filterRight(Ar,B1,B2,C,D1r,D2r,mu1,opts);

% hinf test
P1 = ss(Ar,B1,C,D1r,-1);
P2 = ss(Ar,B2,C,D2r,-1);
Fr  = ss(Af,Bf,Cf,Df,-1);
[hinfnorm(P1),hinfnorm(P2),hinfnorm(Fr),hinfnorm(P1*Fr-P2)]

[Af,Bf,Cf,Df] = filterRight(Ar,B1,B2,C,D1r,D2r,mu1,2);
Fr  = ss(Af,Bf,Cf,Df,-1);
[hinfnorm(P1),hinfnorm(P2),hinfnorm(Fr),hinfnorm(P1*Fr-P2)]

fprintf('Paused: type any key to continue \n\n ')
pause

%% Left filtering
m   = 1; p1  = m2; p2  = m1;
Al  = rand(n);
Al  = Al/max(abs(eig(Al)))/1.1;  % stable A

B  = rand(n,m);
C1 = rand(p1,n);
C2 = rand(p2,n);
D1l = zeros(p1,m);
D2l = zeros(p2,m);

% H = [H1 H2] = [A  | B ]
%               [---|-- ]
%               [C1 | D1] 
%               [C2 | D2]
%|F*H1 - H2|_infty^2 < mu

H.A = Al;H.C1 = C1;  H.C2 = C2;
H.B = B; H.D1 = D1l; H.D2 = D2l;

fprintf('Left filtering with standard and extended Hinf LMIs \n ')
mu2 = 1;
[Af,Bf,Cf,Df] = filterLeft(Al,B,C1,C2,D1l,D2l,mu2,1);

% hinf test

H1 = ss(Al,B,C1,D1l,-1);
H2 = ss(Al,B,C2,D2l,-1);
Fl = ss(Af,Bf,Cf,Df,-1);
[hinfnorm(H1),hinfnorm(H2),hinfnorm(Fl),hinfnorm(Fl*H1-H2)]

[Af,Bf,Cf,Df] = filterLeft(Al,B,C1,C2,D1l,D2l,mu2,2);
Fl = ss(Af,Bf,Cf,Df,-1);
[hinfnorm(H1),hinfnorm(H2),hinfnorm(Fl),hinfnorm(Fl*H1-H2)]

fprintf('Paused: type any key to continue \n\n ')
pause
%% Joint filtering
fprintf('Joint filtering with standard and extended Hinf LMIs \n ')

[Af,Bf,Cf,Df,mu1,mu2] = filterJoint(P,H,1);
Fjoint  = ss(Af,Bf,Cf,Df,-1);

% hinf test
[hinfnorm(Fl*H1-H2),hinfnorm(Fjoint*H1-H2),mu2^(0.5)]
[hinfnorm(P1*Fr-P2),hinfnorm(P1*Fjoint-P2),mu1^(0.5)]

fprintf('Hinf cost with standard Hinf LMI: %6.4f \n', mu1 + mu2);


[Af,Bf,Cf,Df,mu1,mu2] = filterJoint(P,H,2);
Fjoint  = ss(Af,Bf,Cf,Df,-1);
% hinf test
[hinfnorm(Fl*H1-H2),hinfnorm(Fjoint*H1-H2),mu2^(0.5)]
[hinfnorm(P1*Fr-P2),hinfnorm(P1*Fjoint-P2),mu1^(0.5)]

fprintf('Hinf cost with extended Hinf LMI: %6.4f \n', mu1 + mu2);
