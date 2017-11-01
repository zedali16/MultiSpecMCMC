function[Beta_mean,Beta_var,yobs_tmp] = postBeta2(chol_index, Phi_temp, j, yobs,...
                        tau_temp_1, tau_temp_2, Beta_temp_1, Beta_temp_2, xi_temp, nseg)
    global dimen options nbasis nBeta sigmasqalpha 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the mean and variance of normal approximation for coefficient
% of basis fucntions that are the same across segments
%
%   Input:
%       1) chol_index - s index matrix
%       2) Phi_temp - indicate variable that determine which coefficient
%       should be different
%       3) j - the first segment
%       4) yobs - time series observations
%       5) tau_temp_1 - current smoothing parameters for the first segment
%       6) tau_temp_2 - current smoothing parameters for the second segment
%       7) Beta_temp_1 - current coefficients for the first segment
%       8) Beta_temp_2 - current coefficients for the second segment
%       9) xi__temp - current partitions
%       10) nseg - number of observation in the first segment
%   Main Outputs:
%       1) Beta_mean - mean of approximated distribution
%       2) Beta_var - variance of approximated distribution
%       3) yobs_tmp - time series observations in selected segment
%
%   Required programs:  Beta_derive2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %pick right portion of the data
    if j==1
        yobs_tmp = yobs(1:xi_temp(j+1),:);
    elseif j==length(xi_temp)-1
        yobs_tmp = yobs(xi_temp(j-1)+1:end,:);
    else
        yobs_tmp = yobs(xi_temp(j-1)+1:xi_temp(j+1),:);
    end
    
    aa = zeros(2^(dimen^2),1);
    for i=1:2^(dimen^2)
        aa(i)=sum(chol_index(i,:)~=chol_index(Phi_temp,:));
    end
    Phi_inv = find(aa==dimen^2);   
    
    %provide initial values
    x = zeros(chol_index(Phi_inv,:)*repmat(nBeta,dimen^2,1),1);
    
    %optimization process
    [Beta_mean,~,~,~,~,Beta_inv_var] = fminunc(@Beta_derive2, x, options, yobs_tmp, chol_index, Phi_inv,... 
                                                tau_temp_1, tau_temp_2, Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg); 
    %Beta_var = Beta_inv_var\eye(size(Beta_inv_var));
    Beta_var = matpower(Beta_inv_var,-1);
