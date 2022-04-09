function clph2Test
clc;
clear yalmip
% ---------------------------------------------------------------------------- %
%               Mixed traffic control example: Parameters
% ---------------------------------------------------------------------------- %
% The following paramters come from
%   Zheng, Y., Wang, J., & Li, K. (2020). Smoothing traffic flow via control of 
%                          autonomous vehicles. IEEE Internet of Things Journal.

alpha  = 0.6;
beta   = 0.9;
v_max  = 30;
s_st   = 5;
s_go   = 35;
s_star = 20;

% system parameter after linearization  
alpha1 = alpha*v_max/2*pi/(s_go-s_st)*sin(pi*(s_star-s_st)/(s_go-s_st));
alpha2 = alpha+beta;
alpha3 = beta;

% dynamics
P1 = [0 -1;alpha1 -alpha2]; P2 = [0 1;0 alpha3];
A = [P1 zeros(2);P2 P1];
B = blkdiag([0;1],[0;1]);
C = [1 0 0 0;0 0 1 0];

% discretization 
dt = 0.1;
Ad = eye(4)+A*dt; Bd = B*dt; Cd = C;

% ---------------------------------------------------------------------------- %
%               Mixed traffic control example: Design a controller
% ---------------------------------------------------------------------------- %
% H2 synthesis via Closed-loop parameterization
n = 4;          % number of states
p = 2;          % number of outputs
m = 2;          % number of inputs
Q = 1*eye(p);   % performance weights
R = 1*eye(m);


opts.N        = 10;
opts.solver   = 'mosek';
opts.costType = 1;
opts.verbose  = 0;

% SLS
opts.type    = 1;
[Ksls,H2sls,~] = clph2(Ad,Bd,Cd,Q,R,opts);

% IOP
opts.type    = 2;
[Kiop,H2iop,~] = clph2(Ad,Bd,Cd,Q,R,opts);
    
% Mixed I -- mixed sls/iop
opts.type    = 3; 
[KmixI,H2ty3,~] = clph2(Ad,Bd,Cd,Q,R,opts);
     
% Mixed II -- mixed sls/iop 
opts.type    = 4; 
[KmixII,H2ty4,~] = clph2(Ad,Bd,Cd,Q,R,opts);
     
H2optimal = [H2iop, H2sls, H2ty3, H2ty4];
fprintf('\n\n The H2 norm of the closed-loop systems are: \n');
fprintf('               using SLS:       %6.2f\n',H2optimal(1));
fprintf('               using IOP:       %6.2f\n',H2optimal(2));
fprintf('               using Mixed I:   %6.2f\n',H2optimal(3));
fprintf('               using Mixed II:  %6.2f\n',H2optimal(4));

% ---------------------------------------------------------------------------- %
%               Mixed traffic control example: Time domain simulaiton
% ---------------------------------------------------------------------------- %
% parameter
T= 5;
Tn = floor(T/dt);   % number of iterations
dx = zeros(n,Tn);
dy = zeros(p,Tn);
du = zeros(m,Tn);
x0 = [5;0;-2;0];  % initial disturbance
G = ss(Ad,Bd,Cd,[],dt);
Time = (0:Tn)*dt;

% SLS
[x,~,u,~] = dynsim(G,Ksls,dt,dx,dy,du,T,x0);
% figure; subplot(2,1,1); plot(Time,u(1,:)); xlabel('time (s)'); ylabel('Control input'); title('SLS')
%         subplot(2,1,2); plot(Time,u(2,:)); xlabel('time (s)'); ylabel('Control input');
figure; 
subplot(4,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1');  title('SLS')
subplot(4,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');
subplot(4,1,3); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('x_3');
subplot(4,1,4); plot(Time,x(4,:)); xlabel('time (s)'); ylabel('x_4');
   
% IOP
[x,~,u,~] = dynsim(G,Kiop,dt,dx,dy,du,T,x0);
figure; 
subplot(4,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1');  title('IOP')
subplot(4,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');
subplot(4,1,3); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('x_3');
subplot(4,1,4); plot(Time,x(4,:)); xlabel('time (s)'); ylabel('x_4');

% Mix-I parameterization
[x,~,u,~] = dynsim(G,KmixI,dt,dx,dy,du,T,x0);
figure; 
subplot(4,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1');  title('Mix I')
subplot(4,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');
subplot(4,1,3); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('x_3');
subplot(4,1,4); plot(Time,x(4,:)); xlabel('time (s)'); ylabel('x_4');

% Mix-II parameterization
[x,~,u,~] = dynsim(G,KmixII,dt,dx,dy,du,T,x0);
figure; 
subplot(4,1,1); plot(Time,x(1,:)); xlabel('time (s)'); ylabel('x_1');  title('Mix II')
subplot(4,1,2); plot(Time,x(2,:)); xlabel('time (s)'); ylabel('x_2');
subplot(4,1,3); plot(Time,x(3,:)); xlabel('time (s)'); ylabel('x_3');
subplot(4,1,4); plot(Time,x(4,:)); xlabel('time (s)'); ylabel('x_4');


% ---------------------------------------------------------------------------- %
%                                   END
% ---------------------------------------------------------------------------- %
fprintf('\n\nclph2 was successfully tested.\n\n')
end



