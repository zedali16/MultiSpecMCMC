function[PI,nseg_prop,xi_prop,tau_prop,Beta_prop,Phi_prop]=... 
death(chol_index,ts,nexp_curr,nexp_prop,...
	tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,Phi_curr_temp,log_move_curr,log_move_prop)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does the death step in the paper
%
%   Input:
%       1) chol_index - s index matrix
%       2) ts - TxN matrix of time series data
%       3) nexp_curr - current number of segment
%       4) nexp_prop - proposed number of segment
%       5) tau_curr_temp - current smoothing parameters
%       6) xi_curr_temp - current partitions
%       7) nseg_curr_temp - current number of observations in each segment
%       8) Beta_curr_temp - current coefficients
%       9) Phi_curr_temp - which component changed
%       10) log_move_curr - probability: proposed to current
%       11) log_move_prop - probability: current to proposed
%   Main Outputs:
%       1) A - acceptance probability
%       2) nseg_prop - proposed number of observations in each segment
%       3) xi_prop - proposed partitions
%       4) tau_prop - proposed smoothing parameters
%       5) Beta_prop - proposed coefficients
%       6) Phi_prop - proposed indicator variable
%
%   Required programs: postBeta1, Beta_derive_1, whittle_like
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global nobs dimen nBeta nbasis sigmasqalpha tmin tau_up_limit

Beta_prop = zeros(nBeta,dimen^2,nexp_prop);
tau_prop = ones(dimen^2,nexp_prop,1);
nseg_prop = zeros(nexp_prop,1);
xi_prop = zeros(nexp_prop,1);
Phi_prop = zeros(nexp_prop,1);

%***************************************
%Draw a partition to delete
%***************************************
cut_del=unidrnd(nexp_curr-1);

j=0;
for k = 1:nexp_prop
    j = j+1;
    if k==cut_del

        %*************************************************************
		%Calculations related to proposed values	
		%*************************************************************
        xi_prop(k) = xi_curr_temp(j+1);
        tau_prop(:,k) = sqrt(tau_curr_temp(:,j).*tau_curr_temp(:,j+1)); %Combine two taus into one
        nseg_prop(k) = nseg_curr_temp(j) + nseg_curr_temp(j+1); %Combine two segments into one
        Phi_prop(k) = Phi_curr_temp(j+1);
        
        %===================================================================
        %Evaluate the Likelihood at proposed values 
        %===================================================================
        need = sum(Beta_curr_temp(:,:,k)-Beta_curr_temp(:,:,k+1));
        need = need./need;
        need(isnan(need))=0;
        aa = zeros(2^(dimen^2),1);
        for i=1:2^(dimen^2)
            aa(i)=sum(chol_index(i,:)==need);
        end
        Phi_need = find(aa==dimen^2);      
        select = find(chol_index(Phi_need,:)~=0);
        select_inv = find(chol_index(Phi_need,:)==0);
        
        %Compute mean and variances for coefficents of basis functions
        [Beta_mean, Beta_var, yobs_tmp]= postBeta1(chol_index, Phi_need, k, ts, tau_prop(:,k),... 
                                                Beta_curr_temp(:,:,k), xi_prop);
        
        Beta_prop(:,select,k) = reshape(mvnrnd(Beta_mean,0.5*(Beta_var+Beta_var')),nBeta,length(select)); 
        Beta_prop(:,select_inv,k) = Beta_curr_temp(:,select_inv,k);

        %Loglikelihood at proposed values
        [loglike_prop]=whittle_like(yobs_tmp,Beta_prop(:,:,k));
        
        %=============================================================================
		%Evaluate the Proposal Densities at the Proposed values for tau, Phi, and Beta
		%=============================================================================
        pb = reshape(Beta_prop(:,select,k), numel(Beta_prop(:,select,k)),1); 
        %Proposed density for coefficient of basis functions
        log_Beta_prop = - 0.5*(pb-Beta_mean)'*matpower(0.5*(Beta_var+Beta_var'),-1)*(pb-Beta_mean);
                                
        log_seg_prop = -log(nexp_curr-1);  %Proposal for segment choice
        log_Phi_prop = -log(1); %proposal for component choice
        %Calcualte Jacobian
		log_jacobian = -sum(log(2*(sqrt(tau_curr_temp(select,j)) + sqrt(tau_curr_temp(select,j+1))).^2)); 
        
        %log proposal probabililty
        log_proposal_prop = log_Beta_prop + log_seg_prop + log_move_prop + log_Phi_prop;   

        %===========================================================================
		%Evaluate the PRIOR Densities at the Proposed values for tau, Phi, and Beta
		%===========================================================================
        prior_tau = reshape([[repmat(sigmasqalpha,1,dimen^2-dimen*(dimen-1)/2) tau_prop((dimen + dimen*(dimen-1)/2 + 1):end,k)'];...
                        reshape(kron(tau_prop(:,k),ones(nbasis,1)), nbasis,...
                        dimen^2)], nBeta,(dimen^2));  
        prior_tau = reshape(prior_tau(:,select),length(select)*nBeta,1);  
        
        %Prior density for coefficient of basis functions
        log_Beta_prior_prop = -0.5*(pb)'*matpower(diag(prior_tau),-1)*(pb); 
                                            
		%Prior Density of tausq
		log_tau_prior_prop = -length(select)*log(tau_up_limit);
        %Prior Density of Phi
        log_Phi_prior_prop = -log(2^(dimen^2));   
        log_prior_prop = log_tau_prior_prop + log_Beta_prior_prop + log_Phi_prior_prop;
        %*************************************************************
		%Calculations related to current values			
		%*************************************************************
        
		%==================================================================================
		%Evaluate the Likelihood, Proposal and Prior Densities at the Current values
		%==================================================================================
		log_Beta_curr = 0;
		log_tau_prior_curr = 0;
		log_Beta_prior_curr = 0;
		loglike_curr=0;
        for jj=j:j+1
            [Beta_mean, Beta_var, yobs_tmp]= postBeta1(chol_index, Phi_need, jj, ts, tau_curr_temp(:,jj),... 
                                              Beta_curr_temp(:,:,jj), xi_curr_temp);
                                          
             pb = reshape(Beta_curr_temp(:,select,jj), numel(Beta_curr_temp(:,select,jj)),1); 
             %Current density for coefficient of basis functions
             log_Beta_curr = log_Beta_curr -0.5*(pb-Beta_mean)'*matpower(0.5*(Beta_var+Beta_var'),-1)*...
                                                        (pb-Beta_mean);
                                                    
            prior_tau = reshape([[repmat(sigmasqalpha,1,dimen^2-dimen*(dimen-1)/2) tau_curr_temp((dimen + dimen*(dimen-1)/2 + 1):end,jj)'];...
                        reshape(kron(tau_curr_temp(:,jj),ones(nbasis,1)), nbasis,...
                        dimen^2)], nBeta,(dimen^2));  
            prior_tau = reshape(prior_tau(:,select),length(select)*nBeta,1);
            
            %Prior density for coefficient of basis functions at current values                
            log_Beta_prior_curr = log_Beta_prior_curr - 0.5*(pb)'*matpower(diag(prior_tau),-1)*(pb); 
            [log_curr_spec_dens] = whittle_like(yobs_tmp,Beta_curr_temp(:,:,jj)); 
            %Loglikelihood at proposed values
            loglike_curr = loglike_curr + log_curr_spec_dens;
            
            %prior density for smoothing parameters
            log_tau_prior_curr = log_tau_prior_curr - length(select)*log(tau_up_limit);
        end
        
        log_Phi_curr = -log(2^(dimen^2)-1); %proposal for component choice
        log_Phi_prior_curr = -log(2^(dimen^2)); %prior for component choice
        
        %Calculate Log proposal density at current values
        log_proposal_curr = log_move_curr + log_Beta_curr + log_Phi_curr;
        
        %Calculate Priors at Current Vlaues
        log_prior_curr = log_Beta_prior_curr + log_tau_prior_curr + log_Phi_prior_curr;
        j=j+1;
    else
        xi_prop(k) = xi_curr_temp(j);
		tau_prop(:,k) = tau_curr_temp(:,j);
		nseg_prop(k) = nseg_curr_temp(j);
        Beta_prop(:,:,k) = Beta_curr_temp(:,:,j);
        Phi_prop(k) = Phi_curr_temp(j);
    end
end

%=================================================
%Evaluate Target density at proposed values
%=================================================
log_prior_cut_prop=0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop=-log(nobs-(nexp_prop-k+1)*tmin+1);
	else
		log_prior_cut_prop=log_prior_cut_prop-log(nobs-xi_prop(k-1)-(nexp_prop-k+1)*tmin+1);
	end
end
log_target_prop = loglike_prop + log_prior_prop + log_prior_cut_prop;

%==================================================
%Evaluate Target density at current values
%==================================================
log_prior_cut_curr=0;
for k=1:nexp_curr-1
	if k==1
		log_prior_cut_curr=-log(nobs-(nexp_curr-k+1)*tmin+1);
	else
		log_prior_cut_curr=log_prior_cut_curr-log(nobs-xi_curr_temp(k-1)-(nexp_curr-k+1)*tmin+1);
	end
end        
log_target_curr = loglike_curr + log_prior_curr + log_prior_cut_curr;

%*************************************************************
%Calculations acceptance probability	
%*************************************************************
PI = min(1,exp(log_target_prop - log_target_curr + ...
            log_proposal_curr - log_proposal_prop + log_jacobian));        

        
        
        
        
        
        
        
