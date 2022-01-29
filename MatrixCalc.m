function [ Ap,An,As,Ae,Aw ] = MatrixCalc( dx,dy,M,N )

Ap = zeros(1,(M+1)*(N+1));
An = zeros(1,(M+1)*(N+1));
Ae = zeros(1,(M+1)*(N+1));
Aw = zeros(1,(M+1)*(N+1));
As = zeros(1,(M+1)*(N+1));
i = 1;
for j = 1:N+1
    ip = (j-1)*(M+1)+i;
    Ap(ip) = 1.0;
end
i = M+1;
for j = 1:N+1
    ip = (j-1)*(M+1)+i;
    Ap(ip) = 1.0;
end
j = 1;
for i = 2:M
    ip = (j-1)*(M+1)+i;
    Ap(ip) = 1.0;
end
j = N+1;
for i = 2:M
    ip = (j-1)*(M+1)+i;
    Ap(ip) = 1.0;
end
for i = 2:M
    for j = 2:N
        ip = (j-1)*(M+1)+i;
        Aw(ip) = 1.0/(dx^2);
        As(ip) = 1.0/(dy^2);
        Ae(ip) = 1.0/(dx^2);
        An(ip) = 1.0/(dy^2);
        Ap(ip) = (-2.0/(dx^2))+(-2.0/(dy^2));
    end
end

end

