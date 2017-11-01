function[Beta_out, m_out, m, yobs_tmp] = Hamilt2(chol_index, Phi_temp, j, yobs,...
                        tau_temp_1, tau_temp_2, Beta_temp_1, Beta_temp_2, xi_temp, nseg)
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does HMC updates for the coefficents that are the same
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
%       1) Beta_out - vector of coefficients 
%       2) m_out - updated momentum variables
%       3) m - initial momentum variables
%       4) yobs_tmp - time series observations in selected segment
%
%   Required programs:  Gradient2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
    global dimen nbasis sigmasqalpha nBeta M ee
    
    
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
    
    select = chol_index(Phi_inv,:).*(1:dimen^2);
    select = select(select~=0);
    
    ll = nBeta*length(select);
    [gr] = Gradient2(yobs_tmp, chol_index, Phi_inv, tau_temp_1,...
                    tau_temp_2, Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg);
    Beta_old = reshape(Beta_temp_1(:,select),ll,1);
    Beta_out = Beta_old;
    % generate momentum variable
    m = mvnrnd(zeros(ll,1),M*eye(ll))';
    % determine leap number and stepsize
    stepsize = unifrnd(0,2*ee);
    leap = randsample(1:2*ceil(1/ee),1);

    m_out = m + 0.5*gr*stepsize;
    for i=1:leap
        Beta_out = Beta_out + stepsize*(1/M)*eye(ll)*m_out;
        Beta_temp_1(:,select) = reshape(Beta_out,nBeta,length(select));
        Beta_temp_2(:,select) = reshape(Beta_out,nBeta,length(select));
        if i==leap
            [gr] = Gradient2(yobs_tmp, chol_index, Phi_inv, tau_temp_1,...
                        tau_temp_2, Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg);   
            m_out = m_out + 0.5*gr*stepsize;
        else
            [gr] = Gradient2(yobs_tmp, chol_index, Phi_inv, tau_temp_1,...
                        tau_temp_2, Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg);   
            m_out = m_out + 1*gr*stepsize;
        end    
    end    
    m_out = -m_out;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
