function [f,gr,h] = Beta_derive2(x, yobs_tmp, chol_index, Phi_temp, tau_temp_1,...
                    tau_temp_2, Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function used for optimization process for coeeficients are the same
% across segments
%
%   Input:
%       1) x - initial values for coeeficient of basis functions need to be
%       optimized 
%       2) yobs_tmp - time series data within the segment
%       3) chol_index - index matrix
%       4) Phi_temp - which component changed
%       5) tau_temp_1 - smoothing parameters for the first segment
%       6) tau_temp_2 - smoothing parameters for the second segment
%       7) Beta_temp_1 - current coefficients for the first segment
%       8) Beta_temp_2 - current coefficients for the second segment
%       9) sigmasqalpha - smoothing parameters for the constant in real
%       components
%       10) nbasis - number of basis function used
%       11) nseg - number of observation in the first segment
%   Main Outputs:
%       1) f - log posterior probability based on input parameters
%       2) gr - gradients for optimization process
%       3) h - hession matrix for optimization process
%
%   Required programs: lin_basis_func, Beta_derive1
yobs_tmp_1 = yobs_tmp(1:nseg,:);
yobs_tmp_2 = yobs_tmp(nseg+1:end,:);
[f1,grad1,hes1] = Beta_derive1(x, yobs_tmp_1, chol_index, Phi_temp, tau_temp_1, Beta_temp_1, sigmasqalpha, nbasis);
[f2,grad2,hes2] = Beta_derive1(x, yobs_tmp_2, chol_index, Phi_temp, tau_temp_2, Beta_temp_2, sigmasqalpha, nbasis);

f = f1 + f2;
gr = grad1 + grad2;
h = hes1 + hes2;
