function[Beta_mean,Beta_var,yobs_tmp] = postBeta1(chol_index, Phi_temp,...
                                        j, yobs, tau_temp, Beta_temp, xi_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the mean and variance of normal approximation for coefficient
% of basis fucntions that are different across segments
%
%   Input:
%       1) chol_index - s index matrix
%       2) Phi_temp - indicate variable that determine which coefficient
%       should be different
%       3) j - one of two segments
%       4) yobs - time series observations 
%       5) tau_temp - current smoothing parameters
%       6) Beta_temp - current coefficients
%       7) xi__temp - current partitions
%   Main Outputs:
%       1) Beta_mean - mean of approximated distribution
%       2) Beta_var - variance of approximated distribution
%       3) yobs_tmp - time series observations in selected segment
%
%   Required programs:  Beta_derive1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    global nbasis nBeta sigmasqalpha dimen options
    
    %var_inflate_1 is for variance inflation for first real part of Cholesky components; 
    %var_inflate_1 is for variance inflation for rest of components
    
    %pick right portion of the data
    if j>1
            yobs_tmp = yobs((xi_temp(j-1)+1):xi_temp(j),:);
    else
            yobs_tmp = yobs(1:xi_temp(j),:);
    end
    
    %provide initial values
    x = zeros(chol_index(Phi_temp,:)*repmat(nBeta,dimen^2,1),1);
    
    %optimization process
    [Beta_mean,~,~,~,~,Beta_inv_var] = fminunc(@Beta_derive1, x, options, ...
                        yobs_tmp, chol_index, Phi_temp, tau_temp, Beta_temp, sigmasqalpha, nbasis); 
    %Beta_var = Beta_inv_var\eye(size(Beta_inv_var));
    Beta_var = matpower(Beta_inv_var,-1);
    
    