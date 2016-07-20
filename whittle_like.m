function [log_whittle] = whittle_like(yobs_tmp, gamma)
global nbasis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/03/2016 for the MCMC Non-Multi-spectrum analysis
% Calculate local Whittle likelihood
%
%   Input:
%       1) yobs_tmp - time series in the segment
%       2) gamma - coefficient of basis functions
%   Main Outputs:
%       1) log_whitle - log local whittle likelihood
%   Require programs: lin_basis_function

gamma_1 = gamma(:,1:3);
gamma_2 = gamma(1:nbasis,4);

dim = size(yobs_tmp);
n = dim(1);
nfreq = floor(n/2);
tt = (0:nfreq)/(2*nfreq);
yy = fft(yobs_tmp)/sqrt(n);
y = yy(1:(nfreq+1),:);
nf = length(y);

delta_sq = zeros(2,nf);
[xx_r, xx_i]=lin_basis_func(tt);

%theta's
theta_real = xx_r * gamma_1(:,3);
theta_imag = xx_i * gamma_2;
theta = theta_real + sqrt(-1)*theta_imag;

%delta's
delta_sq(1,:) = exp(xx_r * gamma_1(:,1));
delta_sq(2,:) = exp(xx_r * gamma_1(:,2));

if (mod(n,2)==1) %odd n
    log_whittle = -sum(log(delta_sq(1,2:end))' + log(delta_sq(2,2:end))' + conj(y(2:end,1)).*y(2:end,1).*exp(-xx_r(2:end,:)*gamma_1(:,1)) + ...
           conj(y(2:end,2) - theta(2:end).*y(2:end,1)).*(y(2:end,2) - theta(2:end).*y(2:end,1)).*exp(-xx_r(2:end,:)*gamma_1(:,2))) - ...
           0.5*(log(delta_sq(1,1))' + log(delta_sq(2,1))' + conj(y(1,1)).*y(1,1).*exp(-xx_r(1,:)*gamma_1(:,1)) - ...
           conj(y(1,2) - theta(1).*y(1,1)).*(y(1,2) - theta(1).*y(1,1)).*exp(-xx_r(1,:)*gamma_1(:,2)));
else
    log_whittle = -sum(log(delta_sq(1,2:nfreq))' + log(delta_sq(2,2:nfreq))' + conj(y(2:nfreq,1)).*y(2:nfreq,1).*exp(-xx_r(2:nfreq,:)*gamma_1(:,1)) + ...
           conj(y(2:nfreq,2) - theta(2:nfreq).*y(2:nfreq,1)).*(y(2:nfreq,2) - theta(2:nfreq).*y(2:nfreq,1)).*exp(-xx_r(2:nfreq,:)*gamma_1(:,2))) - ...
           0.5*(log(delta_sq(1,1)) + log(delta_sq(2,1)) + conj(y(1,1)).*y(1,1).*exp(-xx_r(1,:)*gamma_1(:,1)) + ...
           conj(y(1,2) - theta(1).*y(1,1)).*(y(1,2) - theta(1).*y(1,1)).*exp(-xx_r(1,:)*gamma_1(:,2))) - ...
           0.5*(log(delta_sq(1,end)) + log(delta_sq(2,end)) + conj(y(end,1)).*y((nfreq+1),1).*exp(-xx_r(end,:)*gamma_1(:,1)) + ...
           conj(y(end,2) - theta(end).*y(end,1)).*(y(end,2) - theta(end).*y(end,1)).*exp(-xx_r(end,:)*gamma_1(:,2)));
end

