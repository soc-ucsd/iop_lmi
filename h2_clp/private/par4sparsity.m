function const = par4sparsity(UYv,UUv,N,m,p,T_struct,R_struct)

%encode sparsity constraints

    const = [];
    z     = sym('z');


    %sparsity of UY
    UYs = sym('UY',[m p*(N+1)]);

    UY_tf = zeros(m,p);
    for t = 1:N+1
        UY_tf = UY_tf + UYs(:,[(t-1)*m+1:t*p])/z^(t-1);
    end

    for i=1:m      %ach1
        for j=1:p
            if(T_struct(i,j)==0)
                fprintf(' sparsity of UY:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/p );
                [num,~] = numden(UY_tf(i,j));
                cc      = coeffs(num,z);
                [A_eq,b_eq]    = equationsToMatrix(cc,vec(UYs));
                A_eqs          = double(A_eq);
                b_eqs          = double(b_eq);
                const    = [const, A_eqs*vec(UYv)== b_eqs];
            end
        end
    end



    %SPARSITY OF UU
    UUs = sym('UU',[m m*(N+1)]);

    UU_tf = zeros(m,m);
    for t = 1:N+1
        UU_tf = UU_tf + UUs(:,[(t-1)*m+1:t*m])/z^(t-1);
    end



    for i=1:m     
        for j=1:m
            if(R_struct(i,j)==0)
                fprintf(' sparsity of UU:  Percentage %6.4f \n', 100*(m*(i-1)+j)/m/m );
                [num,~] = numden(UU_tf(i,j));
                cc      = coeffs(num,z);
                [A_eq,b_eq]    = equationsToMatrix(cc,vec(UUs));
                A_eqs          = double(A_eq);
                b_eqs          = double(b_eq);
                const    = [const, A_eqs*vec(UUv)== b_eqs];
            end
        end
    end

end
