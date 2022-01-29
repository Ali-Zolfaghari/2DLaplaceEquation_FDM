function [ Rh ] = ResidualCalc( U,B,Ap,An,As,Ae,Aw,M,N )

Rh = zeros(1,(M+1)*(N+1));

for i = 2:M
    for j = 2:N
        ip = (j-1)*(M+1)+i;
        Sum = Ap(ip)*U(ip)+Aw(ip)*U(ip-1)+As(ip)*U(ip-M-1)+Ae(ip)*U(ip+1)+An(ip)*U(ip+M+1);
        Rh(ip) = B(ip)-Sum;
    end
end

end

