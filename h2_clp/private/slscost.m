function cost = slscost(CRv,CMv,CNv,CLv,Q,R,B,C,T,opts)
   %H2cost in SLS
   p = size(C,1);
   n = size(CRv,1);
   m = size(R,1);
    
    if  opts.costType == 1  % penalize Y, U, W, Z
        cost = trace((CNv(:,1:p)'*C'+eye(p))*Q*(C*CNv(:,1:p)+eye(p))) + trace(CLv(:,1:m)'*R*CLv(:,1:m));
        cost = cost + trace((B'*CRv(:,1:n)'*C')*Q*(C*CRv(:,1:n)*B)) + trace((B'*CMv(:,1:n)'+eye(m))*R*(CMv(:,1:n)*B + eye(m)));
        Qc = C'*Q*C;
        for t = 2:T + 1
            cost = cost + trace(CNv(:,(t-1)*p+1:t*p)'*Qc*CNv(:,(t-1)*p+1:t*p)) + ...
                    trace(CLv(:,(t-1)*m+1:t*m)'*R*CLv(:,(t-1)*m+1:t*m)) + ...
                    trace((B'*CRv(:,(t-1)*n+1:t*n)')*Qc*(CRv(:,(t-1)*n+1:t*n)*B)) + ...
                    trace(B'*CMv(:,(t-1)*n+1:t*n)'*R*CMv(:,(t-1)*n+1:t*n)*B);
        end    
    elseif opts.costType == 2 % only feasibility
        cost = 0;
    end
        
end