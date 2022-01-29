function [ U,Residual,Iteration ] = GSSolver( UII,UI,B,Ap,An,As,Ae,Aw,M,N,Ws,MAXITER,MAXERROR,TYPE )

Iteration = 1;
Residual = 1;

if( TYPE == 1)
    while( Iteration <= MAXITER )
        for i = 2:M
            for j = 2:N
                ip = (j-1)*(M+1)+i;
                Sum = Aw(ip)*UII(ip-1)+As(ip)*UII(ip-M-1)+Ae(ip)*UI(ip+1)+An(ip)*UI(ip+M+1);
                UII(ip) = (B(ip)-Sum)/Ap(ip);
            end
        end
        Residual = sqrt(sum((UII-UI).^2))/((M+1)*(N+1));
        UI = Ws*UII+(1.0-Ws)*UI;
        Iteration = Iteration+1;
    end
    
elseif( TYPE == 2)
    while( Residual >= MAXERROR )
        for i = 2:M
            for j = 2:N
                ip = (j-1)*(M+1)+i;
                Sum = Aw(ip)*UII(ip-1)+As(ip)*UII(ip-M-1)+Ae(ip)*UI(ip+1)+An(ip)*UI(ip+M+1);
                UII(ip) = (B(ip)-Sum)/Ap(ip);
            end
        end
        Residual = sqrt(sum((UII-UI).^2))/((M+1)*(N+1));
        UI = Ws*UII+(1.0-Ws)*UI;
        Iteration = Iteration+1;
    end
    
end

U = UII;
end

