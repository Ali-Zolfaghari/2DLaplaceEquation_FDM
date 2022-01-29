function [ U,Residual,Iteration ] = SIPSolver( UII,UI,B,Ap,An,As,Ae,Aw,M,N,MAXITER,MAXERROR,TYPE )

        LUd = zeros(1,(M+1)*(N+1));
        LUb = zeros(1,(M+1)*(N+1));
        LUc = zeros(1,(M+1)*(N+1));
        LUe = zeros(1,(M+1)*(N+1));
        LUf = zeros(1,(M+1)*(N+1));
        ER = zeros(1,(M+1)*(N+1));
        EL = zeros(1,(M+1)*(N+1));
        P = zeros(1,(M+1)*(N+1));
        UM = zeros(1,(M+1)*(N+1));
        
        for i = 1:M+1
            for j = 1:N+1
                ip = (j-1)*(M+1)+i;
                LUb(ip) = As(ip);
                LUc(ip) = Aw(ip);
            end
        end
        LUd(1) = Ap(1);
        LUe(1) = Ae(1)/LUd(1);
        LUf(1) = An(1)/LUd(1);
        for i = 1:M+1
            for j = 1:N+1
                ip = (j-1)*(M+1)+i;
                if( ip ~= 1 )
                    LUd(ip) = Ap(ip)-LUc(ip)*LUe(ip-1);
                    if( ip > M+1 )
                        LUd(ip) = LUd(ip)-LUb(ip)*LUf(ip-M-1);
                    end
                    LUe(ip) = Ae(ip)/LUd(ip);
                    LUf(ip) = An(ip)/LUd(ip);
                end
            end
        end
        for i = 1:M+1
            for j = 1:N+1
                ip = (j-1)*(M+1)+i;
                if( ip > M+1 )
                    EL(ip) = LUb(ip)*LUe(ip-M-1);
                end
                if( ip > 1 )
                    ER(ip) = LUc(ip)*LUf(ip-1);
                end
            end
        end
Iteration = 1;
Residual = 1;
if( TYPE == 1)
    while( Iteration <= MAXITER )
        for i = 1:M+1
            for j = 1:N+1
                ip = (j-1)*(M+1)+i;
                if( ip+M >(M+1)*(N+1) )
                    P(ip) = ER(ip)*UI(ip-M);
                elseif( ip-M < 1 )
                    P(ip) = EL(ip)*UI(ip+M);
                else
                    P(ip) = EL(ip)*UI(ip+M)+ER(ip)*UI(ip-M);
                end
            end
        end
        UM(1) = (B(1)+P(1))/LUd(1);
        for j = 1:N+1
            for i = 1:M+1
                ip = (j-1)*(M+1)+i;
                if( ip ~= 1 )
                    if( ip-M-1 < 1 )
                        UM(ip) = (B(ip)+P(ip)-LUc(ip)*UM(ip-1))/LUd(ip);
                    else
                        UM(ip) = (B(ip)+P(ip)-LUb(ip)*UM(ip-M-1)-LUc(ip)*UM(ip-1))/LUd(ip);
                    end
                end
            end
        end
        UII((M+1)*(N+1)) = UM((M+1)*(N+1));
        for j = N+1:-1:1
            for i = M+1:-1:1
                ip = (j-1)*(M+1)+i;
                if( ip ~= (M+1)*(N+1) )
                    if( ip+M+1 > (M+1)*(N+1) )
                        UII(ip) = UM(ip)-LUe(ip)*UII(ip+1);
                    else
                        UII(ip) = UM(ip)-LUf(ip)*UII(ip+M+1)-LUe(ip)*UII(ip+1);
                    end
                end
            end
        end
        Residual = sqrt(sum((UII-UI).^2))/((M+1)*(N+1));
        UI = UII;
        Iteration = Iteration+1;
    end
elseif( TYPE == 2)
    while( Residual >= MAXERROR )
        for i = 1:M+1
            for j = 1:N+1
                ip = (j-1)*(M+1)+i;
                if( ip+M >(M+1)*(N+1) )
                    P(ip) = ER(ip)*UI(ip-M);
                elseif( ip-M < 1 )
                    P(ip) = EL(ip)*UI(ip+M);
                else
                    P(ip) = EL(ip)*UI(ip+M)+ER(ip)*UI(ip-M);
                end
            end
        end
        UM(1) = (B(1)+P(1))/LUd(1);
        for j = 1:N+1
            for i = 1:M+1
                ip = (j-1)*(M+1)+i;
                if( ip ~= 1 )
                    if( ip-M-1 < 1 )
                        UM(ip) = (B(ip)+P(ip)-LUc(ip)*UM(ip-1))/LUd(ip);
                    else
                        UM(ip) = (B(ip)+P(ip)-LUb(ip)*UM(ip-M-1)-LUc(ip)*UM(ip-1))/LUd(ip);
                    end
                end
            end
        end
        UII((M+1)*(N+1)) = UM((M+1)*(N+1));
        for j = N+1:-1:1
            for i = M+1:-1:1
                ip = (j-1)*(M+1)+i;
                if( ip ~= (M+1)*(N+1) )
                    if( ip+M+1 > (M+1)*(N+1) )
                        UII(ip) = UM(ip)-LUe(ip)*UII(ip+1);
                    else
                        UII(ip) = UM(ip)-LUf(ip)*UII(ip+M+1)-LUe(ip)*UII(ip+1);
                    end
                end
            end
        end
        Residual = sqrt(sum((UII-UI).^2))/((M+1)*(N+1));
        UI = UII;
        Iteration = Iteration+1;
    end
end

U = UII;
end

