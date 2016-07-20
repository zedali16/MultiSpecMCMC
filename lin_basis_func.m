function [xx_r,xx_i]= lin_basis_func(freq_hat)% freq are the frequencies

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/02/2016 for the MCMC Non-Multi-spectrum analysis
% Produces linear basis functions
%
%   Input:
%       1) xx_r - linear basis function for real Cholesky components
%       2) ts - linear basis function for imaginary Cholesky components
%   Main Outputs:
%       1) freq_hat - frequencies used
%

    global nbasis ngamma
    nfreq_hat=length(freq_hat);
    xx_r = ones((nfreq_hat),ngamma);
    xx_i = ones((nfreq_hat),nbasis);
    for j=2:ngamma
        xx_r(:,j) = sqrt(2)* cos(2*pi*(j-1)*freq_hat);
    end
    for j =1:nbasis
        xx_i(:,j) = sqrt(2)*sin(2*pi*j*freq_hat);
    end