function cost = par4cost(XYv,UYv,XUv,UUv,Q,R,N,C,opts)
    %%H2cost
    cost = 0;
    [p,~] = size(Q);
    [m,~] = size(R);
    n = size(C,2);
    
    if opts.costType == 1 % penalize Y, U, W, Z
       
        % step 1
        cost = cost + trace((C*XYv(:,1:p)+eye(p))'*Q*(C*XYv(:,1:p)+eye(p)))+trace(UYv(:,1:m)'*R*UYv(:,1:1*m)) + ...
               trace(XUv(:,1:p)'*C'*Q*C*XUv(:,1:p)) + trace(UUv(:,1:m)'*R*UUv(:,1:m));    
        % the rest of steps  
        Qc = C'*Q*C;
        for t = 2: N+1
            cost = cost + trace(XYv(:,(t-1)*p+1:t*p)'*Qc*XYv(:,(t-1)*p+1:t*p)) + trace(UYv(:,(t-1)*m+1:t*m)'*R*UYv(:,(t-1)*m+1:t*m)) + ...
                trace(XUv(:,(t-1)*p+1:t*p)'*Qc*XUv(:,(t-1)*p+1:t*p)) + trace(UUv(:,(t-1)*m+1:t*m)'*R*UUv(:,(t-1)*m+1:t*m));
        end   
    elseif opts.costType == 2  % only feasibility
        cost = 0;
    end
end
