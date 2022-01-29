function [ R2h ] = RestrictCalc( Rh,M,N )

R2h = zeros(1,((M/2)+1)*((N/2)+1));

[ DltIndx ] = DeleteIndexCalc( M,N );

j = 1;
for i = 1:(M+1)*(N+1)
    if( isempty(find(DltIndx == i, 1)) == 1 )
        R2h(j) = Rh(i);
        j = j+1;
    end
end

end

