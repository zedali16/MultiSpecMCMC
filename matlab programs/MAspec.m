function [spect]=MAspec(beta,sigma,freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take bivariate MA coefficient matrix, output the spectral
% beta: 2x(2xlag) coefficient matrix
% sigma: 2x2 covariance matrix
% freq: freqency where spectral calcuated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dim = size(beta);
len = dim(2);

phispect = zeros(2,2,length(freq));
spect = zeros(2,2,length(freq));

for k=1:length(freq)
    phispect(:,:,k) = eye(2);
    for j=1:(len/2)
        if j==1
            bigmat = beta(:,1:2).* exp(-2*pi*sqrt(-1)*freq(k));
        else
            bigmat = beta(:,(2*j-1):2*j) .* exp(-2*j*pi*sqrt(-1)*freq(k));
        end
        phispect(:,:,k) = phispect(:,:,k) + bigmat;
    end
    spect(:,:,k) = phispect(:,:,k)*sigma*conj(phispect(:,:,k));
end    
