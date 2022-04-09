% -----------------------------------------------------------------------
%
% Numerical examples in the Paper
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


%Num = 6:2:14; % number of agents

Num = 2:1:4; % number of agents


n = 2;
m = 1;
p = 1;
Time            = zeros(length(Num),2);
orderController = zeros(length(Num),1);

for idx = 1:length(Num)
    
    yalmip('clear')
    
% -------------------------------------------------------------------------
%
%       Generate system dynamics
%
% -------------------------------------------------------------------------
    N = Num(idx);
    Ac = cell(N,N);   % matrices for A
    Bc = cell(N,1);   % matrices for B
    Cc = cell(N,1);   % matrices for C

    a = 1;
    gamma = 5;

    for i = 1:N
         Ac{i,i} = [1 1; -1 2];
        if i == 1
            Ac{i,i+1} = 1/gamma*exp(-1/a)*eye(2);
        elseif i == N
            Ac{i-1,i} = 1/gamma*exp(-1/a)*eye(2);
        else
            Ac{i-1,i} = 1/gamma*exp(-1/a)*eye(2); 
            Ac{i,i+1} = 1/gamma*exp(-1/a)*eye(2);
        end        
        Cc{i} = [1 0; 0 1];
        Bc{i} = [0;1];
    end

    for i = 1:N
        if i == 1
            Ad = [Ac{i,i}, Ac{i,i+1}, zeros(2, (N-2)*2)];
        elseif i == N
            Ad = [Ad; zeros(2, (N-2)*2), Ac{i-1,i},Ac{i,i}];
        else
            Ad = [Ad; zeros(2, (i-2)*2), Ac{i-1,i},Ac{i,i},Ac{i,i+1},  zeros(2, (N-i-1)*2)];
        end        
    end

    Bd = blkdiag(Bc{:});
    Cd = blkdiag(Cc{:});

% -------------------------------------------------------------------------
%
%       Decentralized control design - LMI via filtering
%
% -------------------------------------------------------------------------

    n = size(Ad,1);
    p = size(Cd,1);
    m = size(Bd,2);

    D = zeros(p,m);

    % Coprime factorization
    poles = rand(n,1) - 0.5;
    K     = place(Ad,Bd,poles);
    F     = place(Ad',Cd',poles);
    L     = -F';

    Ml = ss(Ad + L*Cd, L, Cd, eye(p),-1);  % G = Ml^(-1)*Nl
    Nl = ss(Ad + L*Cd, Bd+L*D, Cd, D,-1);
    Mr = ss(Ad - Bd*K, Bd, -K, eye(m),-1); % G = Nr*Mr^(-1)
    Nr = ss(Ad - Bd*K, Bd,Cd-D*K, D,-1);


    % Decentralized control
    eps = 1;
    opts.decen = 1;              % decentralized control
    opts.agent = N;              % number of agents
    opts.n = size(Ac{1,1},1);    % dimension of each individual agent
    opts.p = size(Cc{1,1},1);
    opts.m = size(Bc{1,1},2);

    opts.cost = 2;               % regularize the state-space realization
    tic
    [K,X,Y]     = filterControllerDec(Ml,Nl,eps,opts);
    Time(idx,1) = toc;


% -------------------------------------------------------------------------
%
%            Centralized H2 optimal control via SLS + FIR
%
% -------------------------------------------------------------------------

    Q = 1*eye(p);   % performance weights
    R = 1*eye(m);

    opts1.N        = 20;       % FIR length N = 20
    opts1.solver   = 'mosek';  % conic solver
    opts1.verbose  = 1;


    opts1.type     = 1;        % SLS
    opts1.costType = 1;        % H2 cost
    tic
    [Ksls,H2sls,~] = clph2(Ad,Bd,Cd,Q,R,opts1);
    Time(idx,2)    = toc;
    try
        orderController(idx) = size(Ksls.A,1);
    catch
        
    end

end

Time
orderController

