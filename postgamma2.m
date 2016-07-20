function[gamma_mean,gamma_var,yobs_tmp] = postgamma2(chol_index, rho_temp, j, yobs,...
                        tau_temp_1, tau_temp_2, gamma_temp_1, gamma_temp_2, xi_temp, nseg)
    global var_inflate_1 var_inflate_2 options nbasis ngamma sigmasqalpha 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/03/2016 for the MCMC Non-Multi-spectrum analysis
% Calculate the mean and variance of normal approximation for coefficient
% of basis fucntions that are the same across segments
%
%   Input:
%       1) chol_index - s index matrix
%       2) rho_temp - indicate variable that determine which coefficient
%       should be different
%       3) j - the first segment
%       4) yobs - time series observations
%       5) tau_temp_1 - current smoothing parameters for the first segment
%       6) tau_temp_2 - current smoothing parameters for the second segment
%       7) gamma_temp_1 - current coefficients for the first segment
%       8) gamma_temp_2 - current coefficients for the second segment
%       9) xi__temp - current partitions
%       10) nseg - number of observation in the first segment
%   Main Outputs:
%       1) gamma_mean - mean of approximated distribution
%       2) gamma_var - variance of approximated distribution
%       3) yobs_tmp - time series observations in selected segment
%
%   Required programs:  gamma_derive2

    %pick right portion of the data
    if j==1
        yobs_tmp = yobs(1:xi_temp(j+1),:);
    elseif j==length(xi_temp)-1
        yobs_tmp = yobs(xi_temp(j-1)+1:end,:);
    else
        yobs_tmp = yobs(xi_temp(j-1)+1:xi_temp(j+1),:);
    end
    
    aa = zeros(16,1);
    for i=1:16
        aa(i)=sum(chol_index(i,:)~=chol_index(rho_temp,:));
    end
    rho_inv = find(aa==4);   
    
    %apply variance inflation
    var1=var_inflate_1*ones(ngamma, 3*ngamma+nbasis);
    var2=[var_inflate_1*zeros(ngamma,ngamma), var_inflate_2*ones(ngamma,2*ngamma+nbasis)];
    var3=var2;
    var4=[var_inflate_1*zeros(nbasis,ngamma), var_inflate_2*ones(nbasis,2*ngamma+nbasis)];
    var=[var1;var2;var3;var4];
    gr_index = (1:(3*ngamma+nbasis)).*[kron(chol_index(rho_inv,1:3),ones(ngamma,1)'), kron(chol_index(rho_inv,4),ones(nbasis,1)')];
    gr_index = gr_index(find(gr_index~=0));
    var = var(gr_index,gr_index);
    
    %provide initial values
    x = zeros(chol_index(rho_inv,:)*[repmat(ngamma,3,1)' nbasis]',1);
    
    %optimization process
    [gamma_mean,~,~,~,~,gamma_inv_var] = fminunc(@gamma_derive2, x, options, yobs_tmp, chol_index, rho_inv,... 
                                                tau_temp_1, tau_temp_2, gamma_temp_1, gamma_temp_2, sigmasqalpha, nbasis, nseg); 
    
    %apply the inflation
%     gamma_var = var.*eye(length(x))/gamma_inv_var;
    gamma_var = var.*(gamma_inv_var\eye(size(gamma_inv_var)));
    