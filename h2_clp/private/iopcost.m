function cost = iopcost(CYv,CUv,CWv,CZv,Q,R,N,opts)
    %%H2cost
    cost = 0;
    [p,~] = size(Q);
    [m,~] = size(R);

    if opts.costType == 1 % penalize Y,U,W,Z
         for t = 1: N+1
            cost = cost + trace(CYv(:,(t-1)*p+1:t*p)'*Q*CYv(:,(t-1)*p+1:t*p))+trace(CUv(:,(t-1)*m+1:t*m)'*R*CUv(:,(t-1)*m+1:t*m));
            cost = cost + trace(CWv(:,(t-1)*p+1:t*p)'*Q*CWv(:,(t-1)*p+1:t*p))+trace(CZv(:,(t-1)*m+1:t*m)'*R*CZv(:,(t-1)*m+1:t*m));
         end
    elseif opts.costType == 2 % feasibility only
        cost = 0;
    end
end