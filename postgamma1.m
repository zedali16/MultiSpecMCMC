function[gamma_mean,gamma_var,yobs_tmp] = postgamma1(chol_index, rho_temp,...
                                        j, yobs, tau_temp, gamma_temp, xi_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/03/2016 for the MCMC Non-Multi-spectrum analysis
% Calculate the mean and variance of normal approximation for coefficient
% of basis fucntions that are different across segments
%
%   Input:
%       1) chol_index - s index matrix
%       2) rho_temp - indicate variable that determine which coefficient
%       should be different
%       3) j - one of two segments
%       4) yobs - time series observations 
%       5) tau_temp - current smoothing parameters
%       6) gamma_temp - current coefficients
%       7) xi__temp - current partitions
%   Main Outputs:
%       1) gamma_mean - mean of approximated distribution
%       2) gamma_var - variance of approximated distribution
%       3) yobs_tmp - time series observations in selected segment
%
%   Required programs:  gamma_derive1
                                
    global var_inflate_1 var_inflate_2 nbasis ngamma sigmasqalpha options
    
    %var_inflate_1 is for variance inflation for first real part of Cholesky components; 
    %var_inflate_1 is for variance inflation for rest of components
    
    %pick right portion of the data
    if j>1
            yobs_tmp = yobs(xi_temp(j-1)+1:xi_temp(j),:);
    else
            yobs_tmp = yobs(1:xi_temp(j),:);
    end
    
    %apply variance inflation
    var1 = var_inflate_1*ones(ngamma, 3*ngamma+nbasis);
    var2 = [var_inflate_1*zeros(ngamma,ngamma), var_inflate_2*ones(ngamma,2*ngamma+nbasis)];
    var3 = var2;
    var4 = [var_inflate_1*zeros(nbasis,ngamma), var_inflate_2*ones(nbasis,2*ngamma+nbasis)];
    var = [var1;var2;var3;var4];
    
    %determine the dimension of the matrix by rho
    gr_index = (1:(3*ngamma+nbasis)).*[kron(chol_index(rho_temp,1:3),ones(ngamma,1)'),...
               kron(chol_index(rho_temp,4),ones(nbasis,1)')];
    gr_index = gr_index(find(gr_index~=0));
    var = var(gr_index,gr_index);
    
    %provide initial values
    x = zeros(chol_index(rho_temp,:)*[repmat(ngamma,3,1)' nbasis]',1);
    
    %optimization process
    [gamma_mean,~,~,~,~,gamma_inv_var] = fminunc(@gamma_derive1, x, options, ...
                        yobs_tmp, chol_index, rho_temp, tau_temp, gamma_temp, sigmasqalpha, nbasis); 
    
    %apply the inflation
%     gamma_var = var.*eye(length(x))*inv(gamma_inv_var);
      gamma_var = var.*(gamma_inv_var\eye(size(gamma_inv_var)));
    