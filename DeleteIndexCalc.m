function [ DltIndx ] = DeleteIndexCalc( M,N )

k = 1;
for i = 2:2:M+1
    for j = 1:N+1
        DltIndx(k) = (j-1)*(M+1)+i;
        k = k+1;
    end
end
for i = 1:M+1
    for j = 2:2:N+1
        DltIndx(k) = (j-1)*(M+1)+i;
        k =k+1;
    end
end
DltIndx = unique(sort(DltIndx));

end

