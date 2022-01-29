function [ Ih ] = ProlongateCalc( U2h,M,N )

Ih = zeros(1,((2*M)+1)*((2*N)+1));

[ DltIndx ] = DeleteIndexCalc( 2*M,2*N );

j = 1;
for i = 1:(2*M+1)*(2*N+1)
    if( isempty(find(DltIndx == i, 1)) == 1 )
        Ih(i) = U2h(j);
        j = j+1;
    end
end

for i = 2:2:2*M+1
    for j = 1:2:2*N+1
        ip = (j-1)*(2*M+1)+i;
        Ih(ip) = 0.5*(Ih(ip-1)+Ih(ip+1));
    end
end

for i = 1:2*M+1
    for j = 2:2:2*N+1
        ip = (j-1)*(2*M+1)+i;
        Ih(ip) = 0.5*(Ih(ip-2*M-1)+Ih(ip+2*M+1));
    end
end

end