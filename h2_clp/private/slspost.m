function [info,K] = slspost(CRv,CMv,CNv,CLv,C,N,m,n,p,info)
% Post processing for SLS

    info.var.R = value(CRv);
    info.var.M = value(CMv);
    info.var.N = value(CNv);
    info.var.L = value(CLv);
    
    temp = C*info.var.N;
    temp(:,1:p) = temp(1:p) + eye(p);
    K1 = Kstatereal(info.var.L,temp,N);  % state space realizaiton K = L(CN+I)^(-1)
    
    K = Kslsreal(info.var.L,info.var.M,info.var.R,info.var.N,N);     % state space realizaiton K = L - MR^{-1}N 

    %  closed-loop responses from optimization
    z = tf('z');
    Rt = zeros(n,n); Mt = zeros(m,n); Nt = zeros(n,p); Lt = zeros(m,p);
    for k = 1:N+1
        Rt = Rt + info.var.R(:,n*(k-1)+1:n*k)*z^(1-k);  % FIR
        Mt = Mt + info.var.M(:,n*(k-1)+1:n*k)*z^(1-k);  % FIR
        Nt = Nt + info.var.N(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
        Lt = Lt + info.var.L(:,p*(k-1)+1:p*k)*z^(1-k);  % FIR
    end
    info.cl.R = Rt;
    info.cl.M = Mt;
    info.cl.N = Nt;
    info.cl.L = Lt;

    info.K1 = K1;
    info.K = K;
    
    
end

