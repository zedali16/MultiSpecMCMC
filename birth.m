function[A,nseg_prop,xi_prop,tau_prop,gamma_prop,rho_prop,z_tmp]=...
birth(chol_index,ts,nexp_curr,nexp_prop,...
	tau_curr_temp,xi_curr_temp,nseg_curr_temp,gamma_curr_temp,rho_curr_temp,...
    log_move_curr,log_move_prop,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/02/2016 for the MCMC Non-Multi-spectrum analysis
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
%       8) gamma_curr_temp - current coefficients
%       9) rho_curr_temp - which component changed
%       10) log_move_curr - probability: proposed to current
%       11) log_move_prop - probability: current to proposed
%   Main Outputs:
%       1) A - acceptance probability
%       2) nseg_prop - proposed number of observations in each segment
%       3) xi_prop - proposed partitions
%       4) tau_prop - proposed smoothing parameters
%       5) gamma_prop - proposed coefficients
%       6) rho_prop - proposed indicator variable
%
%   Required programs: postgamma1, gamma_derive1, whittle_like


global nobs nbasis ngamma sigmasqalpha tmin tau_up_limit

gamma_prop = zeros(ngamma,4,nexp_prop);
tau_prop = ones(4,nexp_prop,1);
nseg_prop = zeros(nexp_prop,1);
xi_prop = zeros(nexp_prop,1);
rho_prop = zeros(nexp_prop-1,1);

%***************************************
%Drawing segment to split
%***************************************
kk = find(nseg_curr_temp>2*tmin); %Number of segments available for splitting
nposs_seg = length(kk);
seg_cut = kk(unidrnd(nposs_seg)); %Drawing segment to split
nposs_cut = nseg_curr_temp(seg_cut)-2*tmin+1; %Drawing new birthed partition


%************************************************
% Proposing new parameters, gamma, rho, tau, xi
%************************************************
for jj=1:nexp_curr
    if jj<seg_cut %nothing updated or proposed here
        xi_prop(jj) = xi_curr_temp(jj);
        tau_prop(:,jj) = tau_curr_temp(:,jj);
        nseg_prop(jj) = nseg_curr_temp(jj);
        gamma_prop(:,:,jj) = gamma_curr_temp(:,:,jj);
        rho_prop(jj) = rho_curr_temp(jj);
    elseif jj==seg_cut  %updating parameters in the selected paritions
        index = unidrnd(nposs_cut);
        if (seg_cut==1)
            xi_prop(seg_cut)=index+tmin-1;
        else
            xi_prop(seg_cut)=xi_curr_temp(jj-1)-1+tmin+index;
        end
        xi_prop(seg_cut+1) = xi_curr_temp(jj);
        z_tmp = Z;
        z_tmp(xi_prop(seg_cut))=1;
        
        %Determine which Cholesky components should change here
        %Sampling index rho from 2 to 15
        rho_prop(seg_cut) = randsample(2:16,1);
        rho_prop(seg_cut+1) = rho_curr_temp(jj);
        
        %Drawing new tausq
        select = find(chol_index(rho_prop(seg_cut),:)~=0);
        select_inv = find(chol_index(rho_prop(seg_cut),:)==0);
        zz= rand(4,1)'.*chol_index(rho_prop(seg_cut),:);
        uu = zz./(1-zz);
        uu(find(uu==0))= 1;
        tau_prop(:,seg_cut)= tau_curr_temp(:,seg_cut).*uu';
        tau_prop(:,seg_cut+1) = tau_curr_temp(:,seg_cut).*(1./uu)';
        
        %Drawing new values for coefficient of basis function for new
        %birthed segments.
        nseg_prop(seg_cut) = index+tmin-1;
        nseg_prop(seg_cut+1) = nseg_curr_temp(jj)-nseg_prop(seg_cut);
        rho_need = rho_prop(seg_cut);
        for k=jj:(jj+1)
            if k==jj
                [gamma_mean_1, gamma_var_1,yobs_tmp_1] = postgamma1(chol_index, rho_need, ...
                k, ts, tau_prop(:,k), gamma_curr_temp(:,:,jj), xi_prop);
                if chol_index(rho_need,4)==1
                    gamma_prop(:,select,k) = reshape([mvnrnd(gamma_mean_1,...
                        0.5*(gamma_var_1+gamma_var_1')), 0],ngamma,length(select));
                    gamma_prop(:,select_inv,k) = gamma_curr_temp(:,select_inv,jj);
                else    
                    gamma_prop(:,select,k) = reshape(mvnrnd(gamma_mean_1,...
                        0.5*(gamma_var_1+gamma_var_1')),ngamma,length(select)); 
                    gamma_prop(:,select_inv,k) = gamma_curr_temp(:,select_inv,jj);
                end    
            else
                [gamma_mean_2, gamma_var_2,yobs_tmp_2] = postgamma1(chol_index, rho_need, ...
                k, ts, tau_prop(:,k), gamma_curr_temp(:,:,jj), xi_prop);
                if chol_index(rho_need,4)==1
                    gamma_prop(:,select,k) = reshape([mvnrnd(gamma_mean_2,...
                        0.5*(gamma_var_2+gamma_var_2')), 0],ngamma,length(select));
                    gamma_prop(:,select_inv,k) = gamma_curr_temp(:,select_inv,jj);
                else    
                    gamma_prop(:,select,k) = reshape(mvnrnd(gamma_mean_2,...
                        0.5*(gamma_var_2+gamma_var_2')),ngamma,length(select)); 
                    gamma_prop(:,select_inv,k) = gamma_curr_temp(:,select_inv,jj);
                end      
            end
        end  
    else  %nothing updated or proposed here
        xi_prop(jj+1) = xi_curr_temp(jj);
		tau_prop(:,jj+1) = tau_curr_temp(:,jj);
		nseg_prop(jj+1) = nseg_curr_temp(jj);
		gamma_prop(:,:,jj+1) = gamma_curr_temp(:,:,jj);
        rho_prop(jj+1) = rho_curr_temp(jj);
    end
end

%Calculating Jacobian
ja = tau_curr_temp(:,seg_cut)./(zz.*(1-zz))';
ja = ja(ja~=Inf);
log_jacobian = sum(log(2*ja)); 

%*************************************************************
%Calculations related to proposed values
%*************************************************************

%================================================================================
%Evaluating the Likelihood, Proposal and Prior Densities at the Proposed values
%================================================================================
log_gamma_prop=0;
log_tau_prior_prop=0;
log_gamma_prior_prop=0;
loglike_prop=0;

for jj=seg_cut:seg_cut+1
    if jj==seg_cut
        gamma_mean = gamma_mean_1;
        gamma_var = gamma_var_1;
        yobs_tmp = yobs_tmp_1;
    else
        gamma_mean = gamma_mean_2;
        gamma_var = gamma_var_2;
        yobs_tmp = yobs_tmp_2;
    end    
    if chol_index(rho_need,4)==1 
        pg = reshape(gamma_prop(:,select,jj), numel(gamma_prop(:,select,jj)),1); 
        %Proposed density for coefficient of basis functions
        log_gamma_prop = log_gamma_prop - 0.5*(pg(1:(end-1))-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*...
                                        (pg(1:(end-1))-gamma_mean);
        prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                    reshape(kron(tau_prop(select,jj),ones(nbasis,1)), nbasis,...
                    length(select))], length(select)*ngamma,1);
        prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);    
        %Prior density for coefficient of basis functions
        log_gamma_prior_prop = log_gamma_prior_prop -...
		0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                                        (pg(1:(end-1))-zeros(length(prior_tau),1)); 
     else
        pg = reshape(gamma_prop(:,select,jj), numel(gamma_prop(:,select,jj)),1); 
        %Proposed density for coefficient of basis functions
        log_gamma_prop = log_gamma_prop - 0.5*(pg-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*(pg-gamma_mean);
        prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                    reshape(kron(tau_prop(select,jj),ones(nbasis,1)), nbasis,...
                    length(select))], length(select)*ngamma,1);
        %Prior density for coefficient of basis functions        
        log_gamma_prior_prop = log_gamma_prior_prop -...
        0.5*(pg-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*(pg-zeros(length(prior_tau),1)); 
     end
     log_tau_prior_prop = log_tau_prior_prop-length(select)*log(tau_up_limit);   %Prior Density of tausq
     [log_prop_spec_dens]=whittle_like(yobs_tmp,gamma_prop(:,:,jj));
     loglike_prop = loglike_prop+log_prop_spec_dens; %Loglikelihood at proposed values
end
        
log_seg_prop=-log(nposs_seg);%Proposal density for segment choice
log_cut_prop=-log(nposs_cut);%Proposal density for partition choice
log_rho_prop=-log(15); %proposal density for component choice
log_prior_rho_prop = -log(16); %prior density for component choice

%Evaluating prior density for cut points at proposed values
log_prior_cut_prop=0;
for k=1:nexp_prop-1
	if k==1
		log_prior_cut_prop=-log(nobs-(nexp_prop-k+1)*tmin+1);
	else
		log_prior_cut_prop=log_prior_cut_prop-log(nobs-xi_prop(k-1)-(nexp_prop-k+1)*tmin+1);
	end
end
%Calculating Log Proposal density at Proposed values
log_proposal_prop=log_gamma_prop+log_seg_prop+log_move_prop+log_cut_prop+log_rho_prop;
%Calculating Log Prior density at Proposed values
log_prior_prop = log_gamma_prior_prop+log_tau_prior_prop+log_prior_cut_prop+log_prior_rho_prop;
%Calculating Target density at Proposed values
log_target_prop=loglike_prop+log_prior_prop;

%*************************************************************
%Calculations related to current values		
%*************************************************************

%=======================================================================================
%Evaluating the Likelihood, Proposal and Prior Densities at the Current values
%=======================================================================================

[gamma_mean, gamma_var, yobs_tmp]=postgamma1(chol_index, rho_need, seg_cut, ts, tau_curr_temp(:,seg_cut),... 
                                                gamma_curr_temp(:,:,seg_cut), xi_curr_temp);
if chol_index(rho_need,4)==1
    pg = reshape(gamma_curr_temp(:,select,seg_cut),numel(gamma_curr_temp(:,select,seg_cut)),1);
    %Current density for coefficient of basis functions
    log_gamma_curr = - 0.5*(pg(1:(end-1))-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*...
                                        (pg(1:(end-1))-gamma_mean);
    prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                    reshape(kron(tau_curr_temp(select,seg_cut),ones(nbasis,1)),...
                    nbasis, length(select))], length(select)*ngamma,1);
    prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);    
    %Prior density for coefficient of basis functions at current values
    log_gamma_prior_curr = -0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                                        (pg(1:(end-1))-zeros(length(prior_tau),1));
else
    pg = reshape(gamma_curr_temp(:,select,seg_cut),numel(gamma_curr_temp(:,select,seg_cut)),1);
    %Current density for coefficient of basis functions
    log_gamma_curr = -0.5*(pg-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*(pg-gamma_mean);
    prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                    reshape(kron(tau_curr_temp(select,seg_cut),ones(nbasis,1)),...
                    nbasis, length(select))], length(select)*ngamma,1);
    %Prior density for coefficient of basis functions at current values            
    log_gamma_prior_curr =  -0.5*(pg-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                                    (pg-zeros(length(prior_tau),1));
end    
log_tau_prior_curr=-length(select)*log(tau_up_limit); %prior density for smoothing parameters
[log_curr_spec_dens] = whittle_like(yobs_tmp,gamma_curr_temp(:,:,seg_cut));
loglike_curr = log_curr_spec_dens; %Loglikelihood at current values

log_rho_curr = -log(1); %proposal for component choice
%Calculating Log Proposal density at current values
log_proposal_curr = log_gamma_curr + log_move_curr + log_rho_curr;

%Evaluating prior density for partition current values
log_prior_cut_curr=0;
for k=1:nexp_curr-1
	if k==1
		log_prior_cut_curr=-log(nobs-(nexp_curr-k+1)*tmin+1);
	else
		log_prior_cut_curr=log_prior_cut_curr-log(nobs-xi_curr_temp(k-1)-(nexp_curr-k+1)*tmin+1);
	end
end
%Calculating Priors at Current Vlaues
log_prior_curr=log_gamma_prior_curr+log_tau_prior_curr+log_prior_cut_curr;
%Evalulating Target densities at current values
log_target_curr=loglike_curr+log_prior_curr;

%*************************************************************
%Calculations acceptance probability	
%*************************************************************
A = min(1,exp(log_target_prop - log_target_curr+...
              log_proposal_curr - log_proposal_prop + log_jacobian));





        
        


        
        
        
        
        
        
        
        
        