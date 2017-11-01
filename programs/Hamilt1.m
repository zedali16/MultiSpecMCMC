function [Beta_out, m_out, m, yobs_tmp]=Hamilt1(chol_index, Phi_temp,...
                                        j, yobs, tau_temp, Beta_temp, xi_temp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does HMC updates for the coefficents that are different
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
%       1) Beta_out - vector of coefficients 
%       2) m_out - updated momentum variables
%       3) m - initial momentum variables
%       4) yobs_tmp - time series observations in selected segment
%
%   Required programs:  Gradient1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    global nbasis sigmasqalpha dimen nBeta M ee

    %pick right portion of the data
    if j>1
        yobs_tmp = yobs((xi_temp(j-1)+1):xi_temp(j),:);
    else
        yobs_tmp = yobs(1:xi_temp(j),:);
    end
    
    select = chol_index(Phi_temp,:).*(1:dimen^2);
    select = select(select~=0);
    
    
    ll = nBeta*length(select);
    [gr]=Gradient1(yobs_tmp, chol_index, Phi_temp, tau_temp,...
                                    Beta_temp, sigmasqalpha, nbasis);
    Beta_old = reshape(Beta_temp(:,select),ll,1); 
    Beta_out = Beta_old;
    
    % generate momentum variable
    m = mvnrnd(zeros(ll,1),M*eye(ll))';
    % determine leap number and step
    stepsize = unifrnd(0,2*ee);
    leap = randsample(1:2*ceil(1/ee),1);
    %leap = randsample(1:(1/stepsize),1);
    
    m_out = m + 0.5*gr*stepsize;
    for i=1:leap
        Beta_out = Beta_out + stepsize*(1/M)*eye(ll)*m_out;
        Beta_temp(:,select) = reshape(Beta_out,nBeta,length(select));
        if i==leap
            [gr] = Gradient1(yobs_tmp, chol_index, Phi_temp, tau_temp,...
                                    Beta_temp, sigmasqalpha, nbasis);
            m_out = m_out + 0.5*gr*stepsize;
        else
            [gr] = Gradient1(yobs_tmp, chol_index, Phi_temp, tau_temp,...
                                    Beta_temp, sigmasqalpha, nbasis);
            m_out = m_out + 1*gr*stepsize;
        end    
    end    
    m_out = -m_out;
  
