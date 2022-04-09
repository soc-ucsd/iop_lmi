function cost = par3cost(YYv,UYv,UXv,YXv,Q,R,N,B,opts)
    %%H2cost
    cost = 0;
    n = size(B);
    p = size(Q,1);
    m = size(R,1);
     if  opts.costType == 1  % penalize Y, U, W, Z
        for t = 1: N+1
            cost = cost + trace(YYv(:,(t-1)*p+1:t*p)'*Q*YYv(:,(t-1)*p+1:t*p))+trace(UYv(:,(t-1)*m+1:t*m)'*R*UYv(:,(t-1)*m+1:t*m));
            if t == 1
                  cost = cost + trace(B'*YXv(:,1:n)'*Q*YXv(:,1:n)*B) + trace((B'*UXv(:,1:n)'+eye(m))*R*(UXv(:,1:n)*B + eye(m)));
            else
                  cost = cost +  trace(B'*YXv(:,(t-1)*n+1:t*n)'*Q*YXv(:,(t-1)*n+1:t*n)*B) + ...
                    trace(B'*UXv(:,(t-1)*n+1:t*n)'*R*UXv(:,(t-1)*n+1:t*n)*B);
            end
        end
     elseif opts.costType == 2  % only consider feasibility
         cost = 0;
     end
end
