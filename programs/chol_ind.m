function[chol_index]=chol_ind(dimen)

chol_index = zeros(2^(dimen^2), dimen^2 );
k=0;
for i=0:dimen^2
    C = nchoosek(1:dimen^2,i);
    dim = size(C);
    if i==0
        k=k+1;
        chol_index(1,:) = 0;
    else  
        for j=1:dim(1)
            k=k+1;
            chol_index(k,C(j,:)) = 1;
        end
    end    
end

end
