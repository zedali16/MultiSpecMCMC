function[PI,nseg_prop,xi_prop,tau_prop,Beta_prop,Phi_prop]=...
            birth(chol_index,ts,nexp_curr,nexp_prop,...
                    tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,Phi_curr_temp,...
                    log_move_curr,log_move_prop)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does the birth step in the paper
%
%   Input:
%       1) chol_index - index matrix
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
%   Required programs: postBeta1, Beta_derive1, whittle_like
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global nobs dimen nbasis nBeta sigmasqalpha tmin tau_up_limit

Beta_prop = zeros(nBeta,dimen^2,nexp_prop);
tau_prop = ones(dimen^2,nexp_prop,1);
nseg_prop = zeros(nexp_prop,1);
xi_prop = zeros(nexp_prop,1);
Phi_prop = zeros(nexp_prop,1);

%***************************************
%Drawing segment to split
%***************************************
kk = find(nseg_curr_temp>2*tmin); %Number of segments available for splitting
nposs_seg = length(kk);
seg_cut = kk(unidrnd(nposs_seg)); %Drawing segment to split
nposs_cut = nseg_curr_temp(seg_cut)-2*tmin+1; %Drawing new birthed partition

%************************************************
% Proposing new parameters, Beta, Phi, tau, xi
%************************************************
for jj=1:nexp_curr
    if jj<seg_cut %nothing updated or proposed here
        xi_prop(jj) = xi_curr_temp(jj);
        tau_prop(:,jj) = tau_curr_temp(:,jj);
        nseg_prop(jj) = nseg_curr_temp(jj);
        Beta_prop(:,:,jj) = Beta_curr_temp(:,:,jj);
        Phi_prop(jj) = Phi_curr_temp(jj);
    elseif jj==seg_cut  %updating parameters in the selected paritions
        index = unidrnd(nposs_cut);
        if (seg_cut==1)
            xi_prop(seg_cut)=index+tmin-1;
        else
            xi_prop(seg_cut)=xi_curr_temp(jj-1)-1+tmin+index;
        end
        xi_prop(seg_cut+1) = xi_curr_temp(jj);

        %Determine which Cholesky components should change here
        Phi_prop(seg_cut) = randsample(2:2^(dimen^2),1);
        Phi_prop(seg_cut+1) = Phi_curr_temp(jj);
        
        %Drawing new tausq
        select = find(chol_index(Phi_prop(seg_cut),:)~=0);
        select_inv = find(chol_index(Phi_prop(seg_cut),:)==0);
        zz = rand(dimen^2,1)'.*chol_index(Phi_prop(seg_cut),:);
        uu = zz./(1-zz);
        uu(find(uu==0))= 1;
        tau_prop(:,seg_cut)= tau_curr_temp(:,seg_cut).*uu';
        tau_prop(:,seg_cut+1) = tau_curr_temp(:,seg_cut).*(1./uu)';

        %Drawing new values for coefficient of basis function for new birthed segments.
        nseg_prop(seg_cut) = index+tmin-1;
        nseg_prop(seg_cut+1) = nseg_curr_temp(jj)-nseg_prop(seg_cut);
        Phi_need = Phi_prop(seg_cut);
        for k=jj:(jj+1)
            if k==jj
                [Beta_mean_1, Beta_var_1,yobs_tmp_1] = postBeta1(chol_index, Phi_need, ...
                k, ts, tau_prop(:,k), Beta_curr_temp(:,:,jj), xi_prop);
                Beta_prop(:,select,k) = reshape(mvnrnd(Beta_mean_1,...
                            0.5*(Beta_var_1+Beta_var_1')),nBeta,length(select)); 
                Beta_prop(:,select_inv,k) = Beta_curr_temp(:,select_inv,jj);
  
            else
                [Beta_mean_2, Beta_var_2,yobs_tmp_2] = postBeta1(chol_index, Phi_need, ...
                k, ts, tau_prop(:,k), Beta_curr_temp(:,:,jj), xi_prop);    
                Beta_prop(:,select,k) = reshape(mvnrnd(Beta_mean_2,...
                            0.5*(Beta_var_2+Beta_var_2')),nBeta,length(select)); 
                Beta_prop(:,select_inv,k) = Beta_curr_temp(:,select_inv,jj);
            end
        end  
    else  %nothing updated or proposed here
        xi_prop(jj+1) = xi_curr_temp(jj);
		tau_prop(:,jj+1) = tau_curr_temp(:,jj);
		nseg_prop(jj+1) = nseg_curr_temp(jj);
		Beta_prop(:,:,jj+1) = Beta_curr_temp(:,:,jj);
        Phi_prop(jj+1) = Phi_curr_temp(jj);
    end
end

%Calculating Jacobian
ja = tau_curr_temp(:,seg_cut)./(zz.*(1-zz))'; ja = ja(ja~=Inf);
log_jacobian = sum(log(2*ja)); 

%*************************************************************
%Calculations related to proposed values
%*************************************************************

%================================================================================
%Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
%================================================================================
log_Beta_prop = 0;
log_tau_prior_prop = 0;
log_Beta_prior_prop = 0;
loglike_prop = 0;

for jj=seg_cut:seg_cut+1
    if jj==seg_cut
        Beta_mean = Beta_mean_1;
        Beta_var = Beta_var_1;
        yobs_tmp = yobs_tmp_1;
    else
        Beta_mean = Beta_mean_2;
        Beta_var = Beta_var_2;
        yobs_tmp = yobs_tmp_2;
    end    
    pb = reshape(Beta_prop(:,select,jj), numel(Beta_prop(:,select,jj)),1); 
    %Proposed density for coefficient of basis functions
    log_Beta_prop = log_Beta_prop - 0.5*(pb-Beta_mean)'*matpower(0.5*(Beta_var+Beta_var'),-1)*(pb-Beta_mean);
    
    %Prior density for coefficient of basis functions 
    prior_tau = reshape([[repmat(sigmasqalpha,1,dimen^2-dimen*(dimen-1)/2) tau_prop((dimen + dimen*(dimen-1)/2 + 1):end,jj)'];...
                     reshape(kron(tau_prop(:,jj),ones(nbasis,1)), nbasis, dimen^2)], nBeta,(dimen^2));  
    prior_tau = reshape(prior_tau(:,select),length(select)*nBeta,1);
         
    log_Beta_prior_prop = log_Beta_prior_prop - 0.5*(pb)'*matpower(diag(prior_tau),-1)*(pb); 

    log_tau_prior_prop = log_tau_prior_prop-length(select)*log(tau_up_limit);   %Prior Density of tausq
    [log_prop_spec_dens] = whittle_like(yobs_tmp,Beta_prop(:,:,jj));
    loglike_prop = loglike_prop + log_prop_spec_dens; %Loglikelihood at proposed values
end
        
log_seg_prop = -log(nposs_seg);%Proposal density for segment choice
log_cut_prop = -log(nposs_cut);%Proposal density for partition choice
log_Phi_prop = -log(2^(dimen^2)-1); %proposal density for component choice
log_prior_Phi_prop = -log(2^(dimen^2)); %prior density for component choice

%Evaluating prior density for cut points at proposed values
log_prior_cut_prop = 0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop=-log(nobs-(nexp_prop-k+1)*tmin+1);
	else
		log_prior_cut_prop=log_prior_cut_prop-log(nobs-xi_prop(k-1)-(nexp_prop-k+1)*tmin+1);
	end
end
%Calculating Log Proposal density at Proposed values
log_proposal_prop = log_Beta_prop + log_seg_prop + log_move_prop + log_cut_prop + log_Phi_prop;
%Calculating Log Prior density at Proposed values
log_prior_prop = log_Beta_prior_prop + log_tau_prior_prop + log_prior_cut_prop + log_prior_Phi_prop;
%Calculating Target density at Proposed values
log_target_prop = loglike_prop + log_prior_prop;

%*************************************************************
%Calculations related to current values		
%*************************************************************

%=======================================================================================
%Evaluating the Likelihood, Proposal and Prior Densities at the Current values
%=======================================================================================

[Beta_mean, Beta_var, yobs_tmp] = postBeta1(chol_index, Phi_need, seg_cut, ts, tau_curr_temp(:,seg_cut),... 
                                                Beta_curr_temp(:,:,seg_cut), xi_curr_temp);
pb = reshape(Beta_curr_temp(:,select,seg_cut),numel(Beta_curr_temp(:,select,seg_cut)),1);

%Current density for coefficient of basis functions
log_Beta_curr = -0.5*(pb-Beta_mean)'*matpower(0.5*(Beta_var+Beta_var'),-1)*(pb-Beta_mean);

%Prior density for coefficient of basis functions at current values 
prior_tau = reshape([[repmat(sigmasqalpha,1,dimen^2-dimen*(dimen-1)/2) tau_curr_temp((dimen + dimen*(dimen-1)/2 + 1):end,seg_cut)'];...
                  reshape(kron(tau_curr_temp(:,seg_cut),ones(nbasis,1)), nbasis, dimen^2)], nBeta,(dimen^2));  
prior_tau = reshape(prior_tau(:,select),length(select)*nBeta,1);
                                         
log_Beta_prior_curr =  -0.5*(pb)'*matpower(diag(prior_tau),-1)*(pb);
                                
log_tau_prior_curr = -length(select)*log(tau_up_limit); %prior density for smoothing parameters
[log_curr_spec_dens] = whittle_like(yobs_tmp,Beta_curr_temp(:,:,seg_cut));
loglike_curr = log_curr_spec_dens; %Loglikelihood at current values

log_Phi_curr = -log(1); %proposal for component choice
%Calculating Log Proposal density at current values
log_proposal_curr = log_Beta_curr + log_move_curr + log_Phi_curr;

%Evaluating prior density for partition current values
log_prior_cut_curr = 0;
for k=1:nexp_curr-1
	if k==1
		log_prior_cut_curr = -log(nobs-(nexp_curr-k+1)*tmin+1);
	else
		log_prior_cut_curr = log_prior_cut_curr-log(nobs-xi_curr_temp(k-1)-(nexp_curr-k+1)*tmin+1);
	end
end
log_prior_Phi_curr = -log(2^(dimen^2)); %prior density for component choice
%Calculating Priors at Current Vlaues
log_prior_curr = log_Beta_prior_curr + log_tau_prior_curr + log_prior_cut_curr + log_prior_Phi_curr;
%Evalulating Target densities at current values
log_target_curr = loglike_curr + log_prior_curr;

%*************************************************************
%Calculations acceptance probability	
%*************************************************************
    PI = min(1,exp(log_target_prop - log_target_curr +...
              log_proposal_curr - log_proposal_prop + log_jacobian));
end    
   




        
        


        
        
        
        
        
        
        
        
        