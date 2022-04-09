function Constraints = par4cons(G,XYv,XUv,UYv,UUv,N,z,A,B)
%achievability constraints for the third parametrization

    [p,m] = size(G);      % system dimension
    n = size(A,1);
    
    XYs = sym('XY',[n p*(N+1)]);
    XUs = sym('XU',[n m*(N+1)]);
    UYs = sym('UY',[m p*(N+1)]);
    UUs = sym('UU',[m m*(N+1)]);
    
    %%Express  as FIR transfer matrices of order N
    XY = zeros(n,p);
    XU = zeros(n,m);
    UY = zeros(m,p);
    UU = zeros(m,m);
    
    for t = 1:N+1
        XY = XY + XYs(:,(t-1)*p+1:t*p)/z^(t-1);
        XU = XU + XUs(:,(t-1)*m+1:t*m)/z^(t-1);
        UY = UY + UYs(:,(t-1)*p+1:t*p)/z^(t-1);
        UU = UU + UUs(:,(t-1)*m+1:t*m)/z^(t-1);
    end


    %%% achievability constraints
    ach1 = (z*eye(n)-A)*XY - B*UY;  %% CAN BE OPTIMIZED similar to sls 
    ach2 = (z*eye(n)-A)*XU - B*UU;  %% CAN BE OPTIMIZED similar to sls
    ach3 = -XY*G + XU - inv(z*eye(n)-A)*B;
    ach4 = -UY * G + UU - eye(m); 
    
    
    Constraints=[];

    for i=1:n      %ach1
        for j=1:p
            fprintf(' ach1:  Percentage %6.4f \n', 100*(p*(i-1)+j)/p/n );
            [num,~] = numden(ach1(i,j));
            cc      = coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(XYs);vec(UYs)]);
            A_eqs          = double(A_eq);
            b_eqs          = double(b_eq);
            Constraints    = [Constraints, A_eqs*[vec(XYv);vec(UYv)]== b_eqs];
        end
    end

    for i=1:n       %ach2
        for j=1:m
            fprintf(' ach2:  Percentage %6.4f \n', 100*(m*(i-1)+j)/n/m );
            [num,~] = numden(ach2(i,j));
            cc      = coeffs(num,z);
            [A_eq,b_eq] = equationsToMatrix(cc,[vec(XUs);vec(UUs)]);
            A_eqs       = double(A_eq);
            b_eqs       = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(XUv);vec(UUv)]== b_eqs];
        end
    end
    
    for i=1:n       %ach3
        for j=1:m
            fprintf(' ach3:  Percentage %6.4f \n', 100*(m*(i-1)+j)/n/m );
            [num,~]=numden(ach3(i,j));
            cc=coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(XUs);vec(XYs)]);
            A_eqs   = double(A_eq);
            b_eqs    = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(XUv);vec(XYv)]== b_eqs];
        end
    end
    
    for i=1:m       %ach4
        for j=1:m
            fprintf(' ach4:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/m );
            [num,~]=numden(ach4(i,j));
            cc=coeffs(num,z);
            [A_eq,b_eq]    = equationsToMatrix(cc,[vec(UYs);vec(UUs)]);
            A_eqs   = double(A_eq);
            b_eqs   = double(b_eq);
            Constraints = [Constraints, A_eqs*[vec(UYv);vec(UUv)]== b_eqs];
        end
    end
end


