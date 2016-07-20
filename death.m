function[A,nseg_prop,xi_prop,tau_prop,gamma_prop,rho_prop,z_tmp]=... 
death(chol_index,ts,nexp_curr,nexp_prop,...
	tau_curr_temp,xi_curr_temp,nseg_curr_temp,gamma_curr_temp,rho_curr_temp,log_move_curr,log_move_prop,Z)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/03/2016 for the MCMC Non-Multi-spectrum analysis
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
%   Required programs: postgamma1, gamma_derive_1, whittle_like

global nobs ngamma nbasis sigmasqalpha tmin tau_up_limit

gamma_prop = zeros(ngamma,4,nexp_prop);
tau_prop = ones(4,nexp_prop,1);
nseg_prop = zeros(nexp_prop,1);
xi_prop = zeros(nexp_prop,1);
rho_prop = zeros(nexp_prop,1);

%***************************************
%Draw a partition to delete
%***************************************
cut_del=unidrnd(nexp_curr-1);

j=0;
for k = 1:nexp_prop
    j = j+1;
    if k==cut_del
        z_tmp = Z;
        z_tmp(xi_curr_temp(k)) = 0;
        %*************************************************************
		%Calculations related to proposed values	
		%*************************************************************
        xi_prop(k) = xi_curr_temp(j+1);
        tau_prop(:,k) = sqrt(tau_curr_temp(:,j).*tau_curr_temp(:,j+1)); %Combine two taus into one
        nseg_prop(k) = nseg_curr_temp(j) + nseg_curr_temp(j+1); %Combine two segments into one
        rho_prop(k) = rho_curr_temp(j+1);
        
        %===================================================================
        %Evaluate the Likelihood at proposed values 
        %===================================================================
   
        %find which components are differend based on previous iterations
        if rho_curr_temp(k+1)==0
            rho_curr_temp(k+1) = rho_curr_temp(k+1) + 1;
        end
        if cut_del==1
            need = chol_index(rho_curr_temp(k),:) + ...
                chol_index(rho_curr_temp(k+1),:)*(Z(xi_curr_temp(k+1))<=Z(xi_curr_temp(k)));
            need = need./need;
            need(isnan(need))=0;
        elseif cut_del==nexp_curr-1
            need = chol_index(rho_curr_temp(k),:) + ...
                chol_index(rho_curr_temp(k-1),:)*(Z(xi_curr_temp(k-1))<=Z(xi_curr_temp(k)));
            need = need./need;
            need(isnan(need))=0;
        else
            need = chol_index(rho_curr_temp(k),:) + ...
                            chol_index(rho_curr_temp(k-1),:)*(Z(xi_curr_temp(k-1))<=Z(xi_curr_temp(k)))+...
                            chol_index(rho_curr_temp(k+1),:)*(Z(xi_curr_temp(k+1))<=Z(xi_curr_temp(k)));
            need = need./need;
            need(isnan(need))=0;
        end    
        aa = zeros(16,1);
        for i=1:16
            aa(i)=sum(chol_index(i,:)==need);
        end
        rho_need = find(aa==4);      
        cc_1 = find(chol_index(rho_need,:)~=0);
        cc_2 = find(chol_index(rho_need,:)==0);
        
        %Compute mean and variances for coefficents of basis functions
        [gamma_mean, gamma_var, yobs_tmp]= postgamma1(chol_index, rho_need, k, ts, tau_prop(:,k),... 
                                                gamma_curr_temp(:,:,k), xi_prop);
        if chol_index(rho_need,4)==1
            gamma_prop(:,cc_1,k) = reshape([mvnrnd(gamma_mean,0.5*(gamma_var+gamma_var')), 0],ngamma,length(cc_1));
            gamma_prop(:,cc_2,k) = gamma_curr_temp(:,cc_2,k);
        else    
            gamma_prop(:,cc_1,k) = reshape(mvnrnd(gamma_mean,0.5*(gamma_var+gamma_var')),ngamma,length(cc_1)); 
            gamma_prop(:,cc_2,k) = gamma_curr_temp(:,cc_2,k);
        end  
        
        %Loglikelihood at proposed values
        [loglike_prop]=whittle_like(yobs_tmp,gamma_prop(:,:,k));
        
        %=============================================================================
		%Evaluate the Proposal Densities at the Proposed values for tau, rho, and gamma
		%=============================================================================
        select = chol_index(rho_need,:).*[1 2 3 4];
        select = select(select~=0);
        if chol_index(rho_need,4)==1 
            pg = reshape(gamma_prop(:,select,k), numel(gamma_prop(:,select,k)),1); 
            %Proposed density for coefficient of basis functions
            log_gamma_prop = - 0.5*(pg(1:(end-1))-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*...
                                    (pg(1:(end-1))-gamma_mean);
        else
            pg = reshape(gamma_prop(:,select,k), numel(gamma_prop(:,select,k)),1); 
            %Proposed density for coefficient of basis functions
            log_gamma_prop = - 0.5*(pg-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*...
                                    (pg-gamma_mean);
        end
        log_seg_prop=-log(nexp_curr-1);  %Proposal for segment choice
        log_rho_prop=-log(1); %proposal for component choice
        %Calcualte Jacobian
		log_jacobian = -sum(log(2*(sqrt(tau_curr_temp(:,j))+sqrt(tau_curr_temp(j+1))).^2)); 
        
        %log proposal probabililty
        log_proposal_prop = log_gamma_prop + log_seg_prop + log_move_prop + log_rho_prop;   

        %========================================================================
		%Evaluate the PRIOR Densities at the Proposed values for tau, rho, and
		%gamma
		%========================================================================
        if chol_index(rho_need,4)==1 
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                            reshape(kron(tau_prop(select,k),ones(nbasis,1)), nbasis,...
                            length(select))], length(select)*ngamma,1);
            prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);    
            %Prior density for coefficient of basis functions
            log_gamma_prior_prop = -0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                                                (pg(1:(end-1))-zeros(length(prior_tau),1)); 
        else
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                            reshape(kron(tau_prop(select,k),ones(nbasis,1)), nbasis,...
                            length(select))], length(select)*ngamma,1);
            %Prior density for coefficient of basis functions
            log_gamma_prior_prop = -0.5*(pg-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                                                (pg-zeros(length(prior_tau),1)); 
        end        
		%Prior Density of tausq
		log_tau_prior_prop=-length(select)*log(tau_up_limit);
        %rho
        %log_rho_prior_prop=-log(16);   
        log_prior_prop =log_tau_prior_prop+log_gamma_prior_prop;
        %*************************************************************
		%Calculations related to current values			
		%*************************************************************
        
		%==================================================================================
		%Evaluate the Likelihood, Proposal and Prior Densities at the Current values
		%==================================================================================
		log_gamma_curr=0;
		log_tau_prior_curr=0;
		log_gamma_prior_curr=0;
		loglike_curr=0;
        for jj=j:j+1
            [gamma_mean, gamma_var, yobs_tmp]= postgamma1(chol_index, rho_need, jj, ts, tau_curr_temp(:,jj),... 
                                              gamma_curr_temp(:,:,jj), xi_curr_temp);
            if chol_index(rho_need,4)==1
                pg = reshape(gamma_curr_temp(:,select,jj), numel(gamma_curr_temp(:,select,jj)),1); 
                %Current density for coefficient of basis functions
                log_gamma_curr = log_gamma_curr - 0.5*(pg(1:(end-1))-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*...
                                                (pg(1:(end-1))-gamma_mean);
                prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                                    reshape(kron(tau_curr_temp(select,jj),ones(nbasis,1)), nbasis,...
                                    length(select))], length(select)*ngamma,1);
                prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);   
                %Prior density for coefficient of basis functions at current values
                log_gamma_prior_curr = log_gamma_prior_curr -...
                0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*....
                                            (pg(1:(end-1))-zeros(length(prior_tau),1)); 
            else
                pg = reshape(gamma_curr_temp(:,select,jj), numel(gamma_curr_temp(:,select,jj)),1); 
                %Current density for coefficient of basis functions
                log_gamma_curr = log_gamma_curr -0.5*(pg-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*...
                                                        (pg-gamma_mean);
                prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                                    reshape(kron(tau_curr_temp(select,jj),ones(nbasis,1)), nbasis,...
                                    length(select))], length(select)*ngamma,1);
                %Prior density for coefficient of basis functions at current values                
                log_gamma_prior_curr = log_gamma_prior_curr -...
                0.5*(pg-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                                        (pg-zeros(length(prior_tau),1)); 
            end
            [log_curr_spec_dens]=whittle_like(yobs_tmp,gamma_curr_temp(:,:,jj)); 
            %Loglikelihood at proposed values
            loglike_curr=loglike_curr + log_curr_spec_dens;
            
            %prior density for smoothing parameters
            log_tau_prior_curr=log_tau_prior_curr-length(select)*log(tau_up_limit);
        end
        
        log_rho_curr=-log(15); %proposal for component choice
        log_rho_prior_curr = -log(16); %prior for component choice
        
        %Calculate Log proposal density at current values
        log_proposal_curr=log_move_curr + log_gamma_curr + log_rho_curr;
        
        %Calculate Priors at Current Vlaues
        log_prior_curr=log_gamma_prior_curr + log_tau_prior_curr + log_rho_prior_curr;
        j=j+1;
    else
        xi_prop(k)=xi_curr_temp(j);
		tau_prop(:,k)=tau_curr_temp(:,j);
		nseg_prop(k)=nseg_curr_temp(j);
        gamma_prop(:,:,k)=gamma_curr_temp(:,:,j);
        rho_prop(k) = rho_curr_temp(j);
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
log_target_prop=loglike_prop + log_prior_prop + log_prior_cut_prop;

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
log_target_curr=loglike_curr+log_prior_curr+log_prior_cut_curr;

%*************************************************************
%Calculations acceptance probability	
%*************************************************************
A = min(1,exp(log_target_prop - log_target_curr + ...
            log_proposal_curr - log_proposal_prop + log_jacobian));        

        
        
        
        
        
        
        