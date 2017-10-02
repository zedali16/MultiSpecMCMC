function [result] = conv_diag(nobs, nBeta, nseg,nexp_curr, Beta_curr, tausq_curr, dimen)

%   conv_diag Calculate convergence diagnostics for beta and tau sq params
%   Get max eigenvalue of beta and tausq matrices across time and cov

idx = cumsum(nseg);
tausq_diag = zeros(nobs,dimen^2);
Beta_diag = zeros(nobs,nBeta,dimen^2);

for i=1:nexp_curr
    if i==1
        for k=1:dimen^2
            tausq_diag(1:idx(i),k) = repmat(tausq_curr(k,i),nseg(i),1);
            
            for j=1:nBeta
                Beta_diag(1:idx(i),j,k)=repmat(Beta_curr(j,k),nseg(i),1); 
            end    
        end   
    else
        for k=1:dimen^2
            tausq_diag((idx(i-1)+1):idx(i),k) = repmat(tausq_curr(k,i),nseg(i),1);
            
            for j=1:nBeta
                Beta_diag((idx(i-1)+1):idx(i),j,k)=repmat(Beta_curr(j,k),nseg(i),1); 
            end    
        end     
    end    
        
end
    %calculate eigenval of symmetric transform
    result = zeros(1,dimen^2*(nBeta+1));
    for k=1:dimen^2
        result(k) = real(max(eig(tausq_diag(:,k)*tausq_diag(:,k)')));
    end
    k=1;
    for i=1:dimen^2
        for j=1:nBeta
            result(k+dimen^2) = real(max(eig(Beta_diag(:,j,i)*Beta_diag(:,j,i)')));
            k=k+1;
        end
    end    
end