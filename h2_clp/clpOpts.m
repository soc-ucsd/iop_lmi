function options = clpOpts
% default parameters in clph2.m

    options.N       = 8;    % default FIR length
    options.solver  = 'mosek';
    options.verbose = 1;    % display solver information
    options.polish  = 0;    % remove tiny numbers in controller
    options.eps     = 1e-8; % rounding precision 1e-8;
    
   
    options.type    = 2;    % which parameterization is used
                            % 1: SLP; 2: IOP; 3: mixed I, 4: mixed II
   
    options.spa     = 0;    % No sparsity on controller structure
    options.S       = [];   % no controller structure
    options.T       = [];   % T, and R come from the idea of sparsity inviarance
    options.R       = [];
    options.stable  = 1;    % open-loop stable or not (0/1); for better numerical robustness, 
                            % better to be open-loop stable                           
    options.gamma   = [];   % no robust stability constraint
    options.bisec   = 0;    % no golden section iteration
    
    options.costType = 1;   % which cost function is used
                            % 1: penalize all four outputs, 2: no penality,
                            % just find a initial stabilizing controller
end
