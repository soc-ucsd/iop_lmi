% -----------------------------------------------------------------------
%
% Example 1 in the Paper
% Title:   Convex Parameterization of Stabilizing Controllers and 
%                                   its LMI-based Computation via Filtering
% Authors: Mauricio C de Oliveira, Yang Zheng
%
%-------------------------------------------------------------------------


clc;clear; close all;

addpath('..');
addpath('../utils')
addpath('../h2_clp')
addpath('../h2_clp/utils')

% load model data 
load chain-data_ex3

yalmip('clear')

% -------------------------------------------------------------------------
%
%       Method 1: decentralized control via LMI + filtering
%
% -------------------------------------------------------------------------

% coprime factorization

%poles = rand(n,1) - 0.5;
K = place(Ad,Bd,poles);
F = place(Ad',Cd',poles);
L = -F';

Ml = ss(Ad + L*Cd, L, Cd, eye(p),-1);  % G = Ml^(-1)*Nl
Nl = ss(Ad + L*Cd, Bd+L*D, Cd, D,-1);


% Decentralized control design via filtering
eps        = 1;
opts.decen = 1;
opts.agent = N;
opts.n = size(Ac{1,1},1);
opts.p = size(Cc{1,1},1);
opts.m = size(Bc{1,1},2);

opts.cost = 2;    % regularize the controller parameters
tic
[K,X,Y] = filterControllerDec(Ml,Nl,eps,opts);
t_lmi = toc;

%  Closed-loop system
G = ss(Ad,Bd,Cd,D,-1);
H = feedback(G,-K);

fprintf('maximal open-loop poles  : %6.2f \n',max(abs(pole(G))))
fprintf('maximal closed-loop poles: %6.2f \n',max(abs(pole(H))))

% Time domain simulation
dt = 1;
T  = 25*dt;
Tn = floor(T/dt);   % number of iterations
dx = zeros(n,Tn);
dy = zeros(p,Tn);
du = zeros(m,Tn);

x0 = [0;1;0;1;0;1]; % inital state

Time = (0:Tn)*dt;

[x,y,u,~] = dynsim(G,K,dt,dx,dy,du,T,x0);

fontsize = 10;
colorSet = {'r','m','b'};
h = [];
figure
for i = 1:N
    h{i} = plot(Time,y(i,:),colorSet{i},'linewidth',1.2); hold on
end
xlim([0 25])
xlabel('Time step $t$','Interpreter','latex','FontSize',fontsize); 
ylabel('Output measurement','Interpreter','latex','FontSize',fontsize);

h1 = legend('Subsystem 1','Subsystem 2','Subsystem 3','Location','NorthEast');
set(h1,'box','off','Interpreter','latex','FontSize',fontsize)

set(gca,'TickLabelInterpreter','latex');
set(gcf,'Position',[250 150 250 250]);
%print(gcf,'lmi-output','-painters','-depsc2','-r600')
  
figure
for i = 1:N
    h{i} = plot(Time,u(i,:),colorSet{i},'linewidth',1.2); hold on
end
xlim([0 25])
xlabel('Time step $t$','Interpreter','latex','FontSize',fontsize); 
ylabel('Control input','Interpreter','latex','FontSize',fontsize);
h1 = legend('Subsystem 1','Subsystem 2','Subsystem 3','Location','NorthEast');
set(h1,'box','off','Interpreter','latex','FontSize',fontsize)

set(gca,'TickLabelInterpreter','latex');
set(gcf,'Position',[250 150 250 250]);
%print(gcf,'lmi-input','-painters','-depsc2','-r600')

% -------------------------------------------------------------------------
%
%       Method 2:  H2 optimal via SLS using FIR truncation
%
% -------------------------------------------------------------------------

Q = 1*eye(p);             % performance weights
R = 1*eye(m);

opts.N        = 40;       % FIT length
opts.solver   = 'mosek';
opts.verbose  = 1;

opts.type     = 1;        % SLS
opts.costType = 1;        % h2 performance
tic
[Ksls,H2sls,~] = clph2(Ad,Bd,Cd,Q,R,opts);
t_sls = toc;

[t_lmi, t_sls]

% Time-domain simulation for SLS
[x,y,u,kesi] = dynsim(G,Ksls,dt,dx,dy,du,T,x0);

figure
for i = 1:N
    h{i} = plot(Time,y(i,:),colorSet{i},'linewidth',1.2); hold on
end
xlim([0 25])
xlabel('Time step $t$','Interpreter','latex','FontSize',fontsize); 
ylabel('Output measurement','Interpreter','latex','FontSize',fontsize);

h1 = legend('Subsystem 1','Subsystem 2','Subsystem 3','Location','NorthEast');
set(h1,'box','off','Interpreter','latex','FontSize',fontsize)

set(gca,'TickLabelInterpreter','latex');
set(gcf,'Position',[250 150 250 250]);
%print(gcf,'sls-output10','-painters','-depsc2','-r600')
