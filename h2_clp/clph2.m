function [K,H2,info] = clph2(A,B,C,Q,R,userOpts)
%
% clph2: Close-loop parameterization of H2 optimal control for discrete-time systems
%
% Input variables
%      (A,B,C):    system dynamics in discrete time
%      Q:    performance weights on output y
%      R:    performance weights on input u
%
% userOpts is a structure and contains the following options
%      N:      Oder of FIR approximation    (default:8)
%      solver: sedumi, sdpt3, csdp or mosek (default)
%      spa:    Distributed control Yes/No   (default: 0)
%      S:      Sparsity pattern for the controller  (default: [])

% The plant dynamics are:
%                           x = Ax_t + Bu_t + v_t
%                           y = Cx_t + w_t
% The orginal problem is as follows
%                      min_{K}    ||Q^{1/2}Y||^2 + ||R^{1/2}U||^2
%                                    + ||Q^{1/2}W||^2 + ||R^{1/2}Z||^2
%                  subject to     K internally stabilizes G
%                                 K \in S
% where Y, U denote the closed-loop response from w to y and u, respectively; 
%       W, Z denote the closed-loop response from v to y and u.
%
% We solve the problem using closed-loop parameterizations, and one of them is
% as follows
%
%              min_{Y,U,W,Z}  ||Q^{1/2}Y||^2 + ||R^{1/2}U||^2
%                                 + ||Q^{1/2}W||^2 + ||R^{1/2}Z||^2
%               subjec to      [I -G][Y W]
%                                    [U Z]  = [I 0]   (1)
%                                 [Y W][-G] = [I]
%                                 [U Z][I]    [0]     (2)
%                              Y,U,W,Z \in FIR(N)     (3)
%                              Y \in R,  U \in T      (4)


% Rely on YALMIP to reformulate the above problem into an SDP, then call
% Mosek/SeDuMi to get a soluton
%
% Authors: Y. Zheng & L. Furieri

% ============================================
%            set uers' options
% ============================================
opts = clpOpts;
if nargin < 5
    error('Robustiop requires three inputs: A, B, C, Q, R; check your inputs.')
elseif nargin > 3
    opts = setUserOpts(opts,userOpts);
end
N     = opts.N;       % FIR horizon
Type  = opts.type;    % which parameterization is used
                      % 1, SLS
                      % 2, IOP
                      % 3, mixed of SLS/IOP
                      % 4, mixed of SLS/IOP

% ========================================================================
%           Dynamics in symbolic form
% ========================================================================
% To do: we need to avoid this symbolic operation
z  = sym('z');
n  = size(A,1);
%G  = C*(z*eye(n)- A)^(-1)*B;  % symbolic form
[p] = size(C,1);              % system dimension
m = size(B,2);

if max(abs(eig(A))) < 1  % open-loop stable systems
    opts.stable = 1;
elseif  max(abs(eig(A))) >= 1 % 
    warning('Closed-loop parameterizations are better suited for open-loop stable plants\n')
    warning('Check closed-loop stability after getting the controller\n')
end

% ========================================================================
%           Controller design using IOP
% ========================================================================
fprintf('==========================================================\n')
fprintf('       H2 synthesis via closed-loop paramterization       \n')
fprintf('==========================================================\n')
fprintf('Number of inputs       : %i\n',m);
fprintf('Number of outputs      : %i\n',p);
fprintf('FIR horizon            : %i\n',N);
fprintf('Parameterization No.   : %i\n',Type);
fprintf('Distributed control    : %i\n',opts.spa);
fprintf('Numerical solver       : %s\n',opts.solver);
fprintf('==========================================================\n');

% ========================================================================
% Start CLP: closed-loop parameterization
% ========================================================================

Step = 1;
switch Type
    case 1    % sls
        % Define FIR variables (3)
        CRv = sdpvar(n,n*(N+1));          % decision variables for R
        CMv = sdpvar(m,n*(N+1));          % decision variables for M
        CNv = sdpvar(n,p*(N+1));          % decision variables for N
        CLv = sdpvar(m,p*(N+1));          % decision variables for L
        
        % achivability constraint (1)-(3)
        fprintf('Step %d: Encoding the achievability constraints ...\n',Step)
        [const] = slscons(A,B,C,CRv,CMv,CNv,CLv,N);
        Step  = Step +1;
        if opts.spa == 1
            fprintf('Step %d: Encoding the sparsity constraints ...\n',Step)
            const_sparsity = slssparsity(C,CNv,CLv,N,m,n,p,opts.T,opts.R);
            const = [const, const_sparsity];
            Step  = Step +1;
        end
        
        % H2 cost
        fprintf('Step %d: Encoding the H2 cost ...\n',Step)
        cost  = slscost(CRv,CMv,CNv,CLv,Q,R,B,C,N,opts);
        Step  = Step +1;
        
    case 2    % iop
        % decision variables
        CYv = sdpvar(p,p*(N+1));          % decision variables for Y
        CUv = sdpvar(m,p*(N+1));          % decision variables for U
        CWv = sdpvar(p,m*(N+1));          % decision variables for W
        CZv = sdpvar(m,m*(N+1));          % decision variables for Z
        
        % achivability constraint (1)-(3)
        fprintf('Step %d: Encoding the achievability constraints ...\n',Step)
        const = iopcons(G,CYv,CUv,CWv,CZv,N,z,opts);
        Step  = Step +1;
        if opts.spa == 1
            fprintf('Step %d: Encoding the sparsity constraints ...\n',Step)
            const_sparsity = iopsparsity(CYv,CUv,N,m,p,opts.T,opts.R);
            const = [const, const_sparsity];
            Step  = Step +1;
        end
        
        % H2 cost
        fprintf('Step %d: Encoding the H2 cost ...\n',Step)
        cost = iopcost(CYv,CUv,CWv,CZv,Q,R,N,opts);
        Step = Step +1;
        
    case 3  % mix of SLS/IOP
        YXv = sdpvar(p,n*(N+1));          % decision variables for \Phi_yx
        YYv = sdpvar(p,p*(N+1));          % decision variables for \Phi_yy
        UXv = sdpvar(m,n*(N+1));          % decision variables for \Phi_ux
        UYv = sdpvar(m,p*(N+1));          % decision variables for \Phi_uy
        
        % achivability constraint (1)-(3)
        fprintf('Step %d: Encoding the achievability constraints ...\n',Step)
        const = par3cons(G,YXv,YYv,UXv,UYv,N,z,A,C);
        Step  = Step +1;
        if opts.spa == 1
            fprintf('Step %d: Encoding the sparsity constraints ...\n',Step)
            const_sparsity = par3sparsity(UYv,YYv,N,m,p,opts.T,opts.R);
            const = [const, const_sparsity];
            Step  = Step +1;
        end        
        % H2 cost
        fprintf('Step %d: Encoding the H2 cost ...\n',Step)
        cost   = par3cost(YYv,UYv,UXv,YXv,Q,R,N,B,opts);
        Step   = Step +1;
        
    case 4 % mix of SLS/IOP
        XYv = sdpvar(n,p*(N+1));          % decision variables for \Phi_xy
        XUv = sdpvar(n,m*(N+1));          % decision variables for \Phi_xu
        UYv = sdpvar(m,p*(N+1));          % decision variables for \Phi_uy
        UUv = sdpvar(m,m*(N+1));          % decision variables for \Phi_uu
        
        % achivability constraint (1)-(3)
        fprintf('Step %d: Encoding the achievability constraints ...\n',Step)
        const = par4cons(G,XYv,XUv,UYv,UUv,N,z,A,B);
        Step  = Step +1;
        if opts.spa == 1
            fprintf('Step %d: Encoding the sparsity constraints ...\n',Step)
            const_sparsity = par4sparsity(UYv,UUv,N,m,p,opts.T,opts.R);
            const = [const, const_sparsity];
            Step  = Step +1;
        end    
        
        % H2 cost
        fprintf('Step %d: Encoding the H2 cost ...\n',Step)
        cost  = par4cost(XYv,UYv,XUv,UUv,Q,R,N,C,opts);
        Step  = Step +1;
        
    otherwise
        error('Unknown parameterizaiton')
end

% Get a solution
fprintf('Step %d: call a solver to get a solution...\n\n',Step);
yalmipOpts = sdpsettings('solver',opts.solver,'verbose',opts.verbose);
sol        =  optimize(const,cost,yalmipOpts);

if sol.problem == 0
    fprintf('Problem is solved, and an h2 stabilizing controller is found\n')
    H2 = sqrt(value(cost));  % value of the H2 norm!
    fprintf('H2 norm of the normal closed-loop system is %6.4f \n\n', H2);
else
    fprintf('\n Numerical issue ... the problem may be not feasible ...\n\n')
    H2 = sqrt(value(cost));
end


% ========================================================================
%                            Post processing
% ========================================================================
info.H2      = H2;  %  H2 norm of the normial system
info.problem = sol.problem;
if sol.problem == 0
    switch Type
        case 1    % sls
            [info,K] = slspost(CRv,CMv,CNv,CLv,C,N,m,n,p,info); 
        case 2    % iop    
            [info,K] = ioppost(CYv,CUv,CWv,CZv,N,m,p,info);  
        case 3
            [info,K] = mixIpost(YXv,YYv,UXv,UYv,N,m,p,info);  
        case 4
            [info,K] = mixIIpost(XYv,UYv,XUv,UUv,N,m,p,info);
        otherwise
            error('Unknown parameterizaiton')
    end
else
    K = [];
end

end
