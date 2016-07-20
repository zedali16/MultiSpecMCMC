function [spect]=ARspec(phi,sigma,freq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take AR coefficient matrix, output the spectral
% right now only take two-dimensional series
% phi: 2x(2xlag) coefficient matrix
% sigma: 2x2 covariance matrix
% order: lag order
% freq: freqency where spectral calcuated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% 
% phi = [0.5 0; 0 0.5];
% sigma = [1 0.9 ; 0.9 1];


dim = size(phi);
len = dim(2);

phispect = zeros(2,2,length(freq));
spect = zeros(2,2,length(freq));

for k=1:length(freq)
    phispect(:,:,k) = eye(2);
    for j=1:(len/2)
        if j==1
            bigmat = phi(:,1:2).* exp(-2*pi*sqrt(-1)*freq(k));
        else
            bigmat = phi(:,(2*j-1):2*j) .* exp(-2*j*pi*sqrt(-1)*freq(k));
        end
        phispect(:,:,k) = phispect(:,:,k) - bigmat;
    end
    spect(:,:,k) = inv(phispect(:,:,k))*sigma*conj(inv(phispect(:,:,k)));
end    


phi = [0.51 0 0 0; 0 -0.3 0 -0.5];
sigma = [1 0.9 ; 0.9 1];
phi2 = [0.3 0 -0.3 0; 0 -0.3 0 0.3];
sigma2 = [1 0.9 0.9 1];