function [gr] = Gradient2(yobs_tmp, chol_index, Phi_temp, tau_temp_1,...
                    tau_temp_2, Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate log gradients for coeeficients selected to be then same
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
%       2) gr - gradients for optimization process
%
%   Required programs: lin_basis_func, Beta_derive1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

yobs_tmp_1 = yobs_tmp(1:nseg,:);
yobs_tmp_2 = yobs_tmp(nseg+1:end,:);
[grad1] = Gradient1(yobs_tmp_1, chol_index, Phi_temp, tau_temp_1, Beta_temp_1, sigmasqalpha, nbasis);
[grad2] = Gradient1(yobs_tmp_2, chol_index, Phi_temp, tau_temp_2, Beta_temp_2, sigmasqalpha, nbasis);

gr = grad1 + grad2;
