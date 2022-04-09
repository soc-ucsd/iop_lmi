function Constraints = par3cons(G,YXv,YYv,UXv,UYv,N,z,A,C)
%achievability constraints for the third parametrization

    [p,m] = size(G);      % system dimension
    n=size(A,1);
    
    YXs = sym('YX',[p n*(N+1)]);
    YYs = sym('YY',[p p*(N+1)]);
    UXs = sym('UX',[m n*(N+1)]);
    UYs = sym('UY',[m p*(N+1)]);
    
    %%Express  as FIR transfer matrices of order N
    YX = zeros(p,n);
    YY = zeros(p,p);
    UX = zeros(m,n);
    UY = zeros(m,p);
    for t = 1:N+1
        YX = YX + YXs(:,[(t-1)*n+1:t*n])/z^(t-1);
        YY = YY + YYs(:,[(t-1)*p+1:t*p])/z^(t-1);
        UX = UX + UXs(:,[(t-1)*n+1:t*n])/z^(t-1);
        UY = UY + UYs(:,[(t-1)*p+1:t*p])/z^(t-1);
    end


    %%%achievability constraints
    ach1 = YX - G*UX - C*inv(z*eye(n)-A);
    ach2 = YY - G*UY-eye(p);

    
    ach3 = YX*(z*eye(n)-A)-YY*C;  %% CAN BE OPTIMIZED similar to sls 
    ach4 = UX*(z*eye(n)-A)-UY*C; %% CAN BE OPTIMIZED similar to sls
    
    
    Constraints=[];

    for i=1:p      %ach1
        for j=1:n
            
            fprintf(' ach1:  Percentage %6.4f \n', 100*(n*(i-1)+j)/p/n );
            [num,~] = numden(ach1(i,j));
            cc      = coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(YXs);vec(UXs)]);
            A_eqs          = double(A_eq);
            b_eqs          = double(b_eq);
            Constraints    = [Constraints, A_eqs*[vec(YXv);vec(UXv)]== b_eqs];
        end
    end

    for(i=1:p)       %ach2
        for(j=1:p)
            fprintf(' ach2:  Percentage %6.4f \n', 100*(p*(i-1)+j)/p/p );
            [num,~] = numden(ach2(i,j));
            cc      = coeffs(num,z);
            [A_eq,b_eq] = equationsToMatrix(cc,[vec(YYs);vec(UYs)]);
            A_eqs       = double(A_eq);
            b_eqs       = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(YYv);vec(UYv)]== b_eqs];
        end
    end
    
    for(i=1:p)       %ach3
        for(j=1:n)
            fprintf(' ach3:  Percentage %6.4f \n', 100*(n*(i-1)+j)/n/p );
            [num,~]=numden(ach3(i,j));
            cc=coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(YXs);vec(YYs)]);
            A_eqs   = double(A_eq);
            b_eqs    = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(YXv);vec(YYv)]== b_eqs];
        end
    end
    
    for(i=1:m)       %ach4
        for(j=1:n)
            fprintf(' ach4:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/n );
            [num,~]=numden(ach4(i,j));
            cc=coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(UXs);vec(UYs)]);
            A_eqs   = double(A_eq);
            b_eqs    = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(UXv);vec(UYv)]== b_eqs];
        end
    end
end


