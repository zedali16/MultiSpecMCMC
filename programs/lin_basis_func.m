function [xx_r,xx_i]= lin_basis_func(freq_hat)% freq are the frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Produces linear basis functions
%
%   Input:
%       1) xx_r - linear basis function for real Cholesky components
%       2) ts - linear basis function for imaginary Cholesky components
%   Main Outputs:
%       1) freq_hat - frequencies used
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global nBeta dimen
    nfreq_hat=length(freq_hat);
    xx_r = ones((nfreq_hat),nBeta);
    xx_i = ones((nfreq_hat),nBeta);
    for j=2:nBeta
        xx_r(:,j) = sqrt(2)* cos(2*pi*(j-1)*freq_hat)/(2*pi*(j-1));
    end
    for j =1:nBeta
        xx_i(:,j) = sqrt(2)*sin(2*pi*j*freq_hat)/(2*pi*j);
    end
