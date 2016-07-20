function[PI,nseg_new,xi_prop,tau_prop,gamma_prop,rho_prop,seg_temp,z_tmp]=...
        within(chol_index,ts,nexp_temp,tau_temp,xi_curr_temp,nseg_curr_temp,gamma_curr_temp,rho_temp,Z)
global nobs nbasis ngamma sigmasqalpha prob_mm1 tmin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written on 07/02/2016 for the MCMC Non-Multi-spectrum analysis
% Does the within-model move in the paper
%
%   Input:
%       1) chol_index - s index matrix
%       2) ts - TxN matrix of time series data
%       3) nexp_temp - number of segments
%       5) tau_temp -  vector of smoothing parameters
%       6) xi_curr_temp - current partitions
%       7) nseg_curr_temp - current number of observations in each segment
%       8) gamma_curr_temp - current vector of coefficients
%       9) rho_temp - which component changed
%       10) Z - together with rho_temp, indicating how components should be 
%            update
%   Main Outputs:
%       1) PI - acceptance probability
%       2) nseg_new - new number of observations in each segment
%       3) xi_prop - proposed partitions
%       4) tau_prop - proposed smoothing parameters
%       5) gamma_prop - proposed coefficients
%       6) rho_prop - proposed indicator variable
%
%   Required programs: postgamma1, postgamma2, gamma_derive1, gamma_derive2, whittle_like

xi_prop=xi_curr_temp;
gamma_prop=gamma_curr_temp;
nseg_new=nseg_curr_temp;
tau_prop = tau_temp;
rho_prop = rho_temp;


if nexp_temp>1
    %*********************************************************
    % If contains more than one segments
    %*********************************************************

    seg_temp = unidrnd(nexp_temp-1);  %Drawing Segment to cut
    u = rand;
    cut_poss_curr = xi_curr_temp(seg_temp);
    nposs_prior = nseg_curr_temp(seg_temp) + nseg_curr_temp(seg_temp+1) - 2*tmin+1;
    
    %Determing if the relocation is a big jump or small jump
    if u<prob_mm1
        if nseg_curr_temp(seg_temp)==tmin && nseg_curr_temp(seg_temp+1)==tmin
            nposs=1; %Number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp)- 1 + new_index;
        elseif nseg_curr_temp(seg_temp)==tmin
			nposs=2; %Number of possible locations for new cutpoint
			new_index=unidrnd(nposs);%Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp)- 1 + new_index;
        elseif nseg_curr_temp(seg_temp+1)==tmin
			nposs=2; %Number of possible locations for new cutpoint
			new_index = unidrnd(nposs); %Drawing index of new cutpoint
			cut_poss_new = xi_curr_temp(seg_temp) + 1 - new_index;
        else
			nposs=3;% Number of possible locations for new cutpoint
			new_index = unidrnd(nposs);%Drawing index of new cutpoint 
			cut_poss_new = xi_curr_temp(seg_temp) - 2 + new_index;
        end
    else
		new_index=unidrnd(nposs_prior);%
		cut_poss_new = sum(nseg_curr_temp(1:seg_temp-1))- 1 + tmin+new_index;
    end
    
    xi_prop(seg_temp)=cut_poss_new;
    if seg_temp>1
        %Number of observations in lower part of new cutpoin
		nseg_new(seg_temp) = xi_prop(seg_temp) - xi_curr_temp(seg_temp-1); 
    else
		nseg_new(seg_temp) = xi_prop(seg_temp);
    end
    %Number of observations in upper part of new cutpoint
	nseg_new(seg_temp+1) = nseg_curr_temp(seg_temp) + nseg_curr_temp(seg_temp+1) - nseg_new(seg_temp);
    
    %***************************
    % moving Z
    %***************************
    z_tmp=Z;
    z_tmp(xi_prop(seg_temp))= z_tmp(xi_curr_temp(seg_temp));
    z_tmp(xi_curr_temp(seg_temp)) = 0; 
    
    %=========================================================================================
    %Evaluating the cut Proposal density for the cut-point at the cureent and proposed values
    %=========================================================================================
    if(abs(cut_poss_new-cut_poss_curr)>1)
        log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior);
        log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior);
    elseif nseg_curr_temp(seg_temp)==tmin && nseg_curr_temp(seg_temp+1)==tmin
        log_prop_cut_prop=0;
        log_prop_cut_curr=0;
    else
        if (nseg_curr_temp(seg_temp)==tmin || nseg_curr_temp(seg_temp+1)==tmin)
           log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior)+log(1/2)+log(prob_mm1);
        else
           log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior)+log(1/3)+log(prob_mm1); 
        end
        if(nseg_new(seg_temp)==tmin || nseg_new(seg_temp+1)==tmin)
           log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior)+log(1/2)+log(prob_mm1);
        else
           log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior)+log(1/3)+log(prob_mm1); 
        end
    end
    %==========================================================
    %find which components were changed before
    %==========================================================
    if nexp_temp==2;
        need = chol_index(rho_temp(seg_temp),:);
        need = need./need;
        need(isnan(need))=0;
    else
        if seg_temp==1
            need = chol_index(rho_temp(seg_temp),:) +...
                   chol_index(rho_temp(seg_temp+1),:)*(Z(xi_curr_temp(seg_temp+1))<=Z(xi_curr_temp(seg_temp)));
            need = need./need;
            need(isnan(need))=0;
        elseif seg_temp==nexp_temp-1
            need = chol_index(rho_temp(seg_temp),:) +... 
                   chol_index(rho_temp(seg_temp-1),:)*(Z(xi_curr_temp(seg_temp-1))<=Z(xi_curr_temp(seg_temp)));
            need = need./need;
            need(isnan(need))=0;
        else
            need = chol_index(rho_temp(seg_temp),:) +...
                   chol_index(rho_temp(seg_temp-1),:)*(Z(xi_curr_temp(seg_temp-1))<=Z(xi_curr_temp(seg_temp)))+...
                   chol_index(rho_temp(seg_temp+1),:)*(Z(xi_curr_temp(seg_temp+1))<=Z(xi_curr_temp(seg_temp)));
            need = need./need;
            need(isnan(need))=0;
        end    
    end
    aa = zeros(16,1);
    for i=1:16
        aa(i)=sum(chol_index(i,:)==need);
    end
    rho_need = find(aa==4);      
    select = chol_index(rho_need,:).*[1 2 3 4];
    select = select(select~=0);
    select_inv = chol_index(rho_need,:).*[1 2 3 4];
    select_inv = find(select_inv==0);
    
    %==========================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the current values
    %==========================================================================
    loglike_curr = 0;
	log_gamma_curr_temp_1=0;
	log_gamma_prior_curr_1=0;    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For coefficient of basis functions that are different across two segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j=seg_temp:seg_temp+1 
        [gamma_mean, gamma_var, yobs_tmp]=...
            postgamma1(chol_index, rho_need, j, ts, tau_temp(:,j), gamma_curr_temp(:,:,j), xi_curr_temp);
        if chol_index(rho_need,4)==1 
            pg = reshape(gamma_curr_temp(:,select,j), numel(gamma_curr_temp(:,select,j)),1); 
            %Current density for coefficient of basis functions
            log_gamma_curr_temp_1 = log_gamma_curr_temp_1 - 0.5*(pg(1:(end-1))-gamma_mean)'*...
                                    inv(0.5*(gamma_var+gamma_var'))*(pg(1:(end-1))-gamma_mean);
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                        reshape(kron(tau_temp(select,j),ones(nbasis,1)), nbasis,...
                        length(select))], length(select)*ngamma,1);
            prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);  
            %Prior density for coefficient of basis functions at current values
            log_gamma_prior_curr_1 = log_gamma_prior_curr_1 -...
                    0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                    (pg(1:(end-1))-zeros(length(prior_tau),1)); 
        else
            pg = reshape(gamma_curr_temp(:,select,j), numel(gamma_curr_temp(:,select,j)),1); 
            %Current density for coefficient of basis functions
            log_gamma_curr_temp_1 = log_gamma_curr_temp_1 -0.5*(pg-gamma_mean)'*...
                                    inv(0.5*(gamma_var+gamma_var'))*(pg-gamma_mean);
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                        reshape(kron(tau_temp(select,j),ones(nbasis,1)), nbasis,...
                        length(select))], length(select)*ngamma,1);
            %Prior density for coefficient of basis functions at current values        
            log_gamma_prior_curr_1 = log_gamma_prior_curr_1 -...
                        0.5*(pg-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                        (pg-zeros(length(prior_tau),1)); 
        end
         %Loglikelihood at current values
         [log_curr_spec_dens]=whittle_like(yobs_tmp,gamma_curr_temp(:,:,j));  
         loglike_curr=loglike_curr+log_curr_spec_dens;
    end
    
    if rho_need ~=16
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For coefficient of basis functions that are the same across two segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [gamma_mean,gamma_var,~] =...
            postgamma2(chol_index, rho_need, seg_temp, ts, tau_temp(:,seg_temp), tau_temp(:,seg_temp+1),...
                       gamma_curr_temp(:,:,seg_temp), gamma_curr_temp(:,:,seg_temp+1), xi_curr_temp, nseg_curr_temp(seg_temp));
        if chol_index(rho_need,4)==0 
            pg = reshape(gamma_curr_temp(:,select_inv,j), numel(gamma_curr_temp(:,select_inv,j)),1); 
            log_gamma_curr_temp_2 =  - 0.5*(pg(1:(end-1))-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*(pg(1:(end-1))-gamma_mean);
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select_inv));...
                        reshape(kron(tau_temp(select_inv,j),ones(nbasis,1)), nbasis,...
                        length(select_inv))], length(select_inv)*ngamma,1);
            prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);    
            log_gamma_prior_curr_2 =  -0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*...
                            inv(diag(prior_tau))*(pg(1:(end-1))-zeros(length(prior_tau),1));
        else
            pg = reshape(gamma_curr_temp(:,select_inv,j), numel(gamma_curr_temp(:,select_inv,j)),1); 
            log_gamma_curr_temp_2 = -0.5*(pg-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*(pg-gamma_mean);
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select_inv));...
                        reshape(kron(tau_temp(select_inv,j),ones(nbasis,1)), nbasis,...
                        length(select_inv))], length(select_inv)*ngamma,1);
            log_gamma_prior_curr_2 = -0.5*(pg-zeros(length(prior_tau),1))'*...
                            inv(diag(prior_tau))*(pg-zeros(length(prior_tau),1)); 
        end
        log_gamma_curr_temp = log_gamma_curr_temp_1 + log_gamma_curr_temp_2;
        log_gamma_prior_curr = log_gamma_prior_curr_1 + log_gamma_prior_curr_2;
    else
        log_gamma_curr_temp = log_gamma_curr_temp_1;
        log_gamma_prior_curr = log_gamma_prior_curr_1;
    end    
    %=================================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the proposed values
    %=================================================================================
    loglike_prop = 0;
	log_gamma_prop_1=0;
	log_gamma_prior_prop_1=0; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For coefficient of basis functions that are different across two segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    yobs_tmp_2=cell(2,1);
    for j=seg_temp:seg_temp+1
        [gamma_mean, gamma_var, yobs_tmp]=postgamma1(chol_index, rho_need, j,...
                                    ts, tau_prop(:,j), gamma_curr_temp(:,:,j), xi_prop);
        if j==seg_temp
            yobs_tmp_2{1}=yobs_tmp;
        else
            yobs_tmp_2{2}=yobs_tmp;
        end
        if chol_index(rho_need,4)==1
            gamma_prop(:,select,j) = reshape([mvnrnd(gamma_mean,0.5*(gamma_var+gamma_var')), 0],ngamma,length(select));
        else    
            gamma_prop(:,select,j) = reshape(mvnrnd(gamma_mean,0.5*(gamma_var+gamma_var')),ngamma,length(select)); 
        end    
        
        if chol_index(rho_need,4)==1 
            pg = reshape(gamma_prop(:,select,j), numel(gamma_prop(:,select,j)),1); 
            %Proposed density for coefficient of basis functions
            log_gamma_prop_1 = log_gamma_prop_1 - 0.5*(pg(1:(end-1))-gamma_mean)'*...
                               inv(0.5*(gamma_var+gamma_var'))*(pg(1:(end-1))-gamma_mean);
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                        reshape(kron(tau_prop(select,j),ones(nbasis,1)), nbasis,...
                        length(select))], length(select)*ngamma,1);
            prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);    
            %Prior density for coefficient of basis functions
            log_gamma_prior_prop_1 = log_gamma_prior_prop_1 -...
                    0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*...
                    inv(diag(prior_tau))*(pg(1:(end-1))-zeros(length(prior_tau),1));
        else
            pg = reshape(gamma_prop(:,select,j), numel(gamma_prop(:,select,j)),1); 
            %Proposed density for coefficient of basis functions
            log_gamma_prop_1 = log_gamma_prop_1 -0.5*(pg-gamma_mean)'*...
                     inv(0.5*(gamma_var+gamma_var'))*(pg-gamma_mean);
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select));...
                        reshape(kron(tau_prop(select,j),ones(nbasis,1)), nbasis,...
                        length(select))], length(select)*ngamma,1);
            %Prior density for coefficient of basis functions        
            log_gamma_prior_prop_1 = log_gamma_prior_prop_1 -...
           0.5*(pg-zeros(length(prior_tau),1))'*...
                    inv(diag(prior_tau))*(pg-zeros(length(prior_tau),1)); 
        end
    end   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For coefficient of basis functions that are the same across two segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if rho_need ~=16
        [gamma_mean,gamma_var,~] = postgamma2(chol_index, rho_need, seg_temp,...
                                   ts, tau_prop(:,seg_temp), tau_prop(:,seg_temp+1),...
                                   gamma_prop(:,:,seg_temp), gamma_prop(:,:,seg_temp+1), xi_prop, nseg_new(seg_temp));
        if chol_index(rho_need,4)==0
            gamma_prop(:,select_inv,seg_temp) = reshape([mvnrnd(gamma_mean,0.5*(gamma_var+gamma_var')), 0],ngamma,...
                                                length(select_inv));
            gamma_prop(:,select_inv,seg_temp+1) = gamma_prop(:,select_inv,seg_temp);
        else    
            gamma_prop(:,select_inv,seg_temp) = reshape(mvnrnd(gamma_mean,0.5*(gamma_var+gamma_var')),ngamma,...
                                                length(select_inv)); 
            gamma_prop(:,select_inv,seg_temp+1) = gamma_prop(:,select_inv,seg_temp);
        end                                            
        if chol_index(rho_need,4)==0 
            pg = reshape(gamma_prop(:,select_inv,seg_temp), numel(gamma_prop(:,select_inv,seg_temp)),1); 
            %Proposed density for coefficient of basis functions
            log_gamma_prop_2 =  - 0.5*(pg(1:(end-1))-gamma_mean)'*...
                                inv(0.5*(gamma_var+gamma_var'))*(pg(1:(end-1))-gamma_mean);
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select_inv));...
                        reshape(kron(tau_prop(select_inv,seg_temp),ones(nbasis,1)), nbasis,...
                        length(select_inv))], length(select_inv)*ngamma,1);
            prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);    
            %Prior density for coefficient of basis functions 
            log_gamma_prior_prop_2 = -0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*...
                                     inv(diag(prior_tau))*(pg(1:(end-1))-zeros(length(prior_tau),1)); 
        else
            pg = reshape(gamma_curr_temp(:,select_inv,seg_temp), numel(gamma_curr_temp(:,select_inv,seg_temp)),1); 
            %Proposed density for coefficient of basis functions
            log_gamma_prop_2 = -0.5*(pg-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*(pg-gamma_mean);
            prior_tau = reshape([repmat(sigmasqalpha,1,length(select_inv));...
                        reshape(kron(tau_prop(select_inv,seg_temp),ones(nbasis,1)), nbasis,...
                        length(select_inv))], length(select_inv)*ngamma,1);
            %Prior density for coefficient of basis functions         
            log_gamma_prior_prop_2 =  -0.5*(pg-zeros(length(prior_tau),1))'*...
                                      inv(diag(prior_tau))*(pg-zeros(length(prior_tau),1)); 
        end
        log_gamma_prop = log_gamma_prop_1 + log_gamma_prop_2;
        log_gamma_prior_prop = log_gamma_prior_prop_1 + log_gamma_prior_prop_2;
    else
        log_gamma_prop = log_gamma_prop_1 ;
        log_gamma_prior_prop = log_gamma_prior_prop_1 ;
    end    
    %Loglikelihood at proposed values
    for j=seg_temp:seg_temp+1
        if j==seg_temp
            [log_curr_spec_dens]=whittle_like(yobs_tmp_2{1},gamma_prop(:,:,j)); 
        else
            [log_curr_spec_dens]=whittle_like(yobs_tmp_2{2},gamma_prop(:,:,j));
        end
        loglike_prop=loglike_prop+log_curr_spec_dens;    
    end

    %proposal density
    log_proposal_curr = log_gamma_curr_temp + log_prop_cut_curr;
    log_proposal_prop = log_gamma_prop + log_prop_cut_prop;
    
    %target density
    log_prior_cut_prop=0;
	log_prior_cut_curr=0;
    for k=1:nexp_temp-1
        if k==1
            log_prior_cut_prop=-log(nobs-(nexp_temp-k+1)*tmin+1);
			log_prior_cut_curr=-log(nobs-(nexp_temp-k+1)*tmin+1);
		else
			log_prior_cut_prop=log_prior_cut_prop-log(nobs-xi_prop(k-1)-(nexp_temp-k+1)*tmin+1);
			log_prior_cut_curr=log_prior_cut_curr-log(nobs-xi_curr_temp(k-1)-(nexp_temp-k+1)*tmin+1);
        end
    end
    log_target_prop=loglike_prop+log_gamma_prior_prop+log_prior_cut_prop;
    log_target_curr=loglike_curr+log_gamma_prior_curr+log_prior_cut_curr;
else    
    
    %*********************************************************
    % If contains only one segment
    %*********************************************************
    
    %%%%%%%%%%%%%%
    %moving Z
    %%%%%%%%%%%%%%
    z_tmp=Z;
    nseg_new=nobs;
    seg_temp=1;
    
    %=================================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the proposed values
    %=================================================================================
    rho_need = 16;
    [gamma_mean, gamma_var, yobs_tmp]=postgamma1(chol_index, rho_need, 1, ts,...
                                        tau_temp, gamma_curr_temp, xi_curr_temp);
    gamma_prop(:,:,1) = reshape([mvnrnd(gamma_mean,0.5*(gamma_var+gamma_var')), 0],ngamma,4);
    
    %Proposed density for coefficient of basis functions
    pg = reshape(gamma_prop(:,:,1), numel(gamma_prop(:,:,1)),1); 
    log_gamma_prop = - 0.5*(pg(1:(end-1))-gamma_mean)'*...
                        inv(0.5*(gamma_var+gamma_var'))*(pg(1:(end-1))-gamma_mean);
    prior_tau = reshape([repmat(sigmasqalpha,1,4);...
                        reshape(kron(tau_prop(:,1),ones(nbasis,1)), nbasis, 4)], 4*ngamma,1);
    prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);  
    %Prior density for coefficient of basis functions
    log_gamma_prior_prop = -0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*inv(diag(prior_tau))*...
                            (pg(1:(end-1))-zeros(length(prior_tau),1));
    %Loglike at proposed values
    [loglike_prop]=whittle_like(yobs_tmp,gamma_prop(:,:,1));
    
    %=================================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the current values
    %=================================================================================
    
    %Proposed density for coefficient of basis functions
    pg = reshape(gamma_curr_temp(:,:,1), numel(gamma_curr_temp(:,:,1)),1); 
    log_gamma_curr_temp = -0.5*(pg(1:(end-1))-gamma_mean)'*inv(0.5*(gamma_var+gamma_var'))*...
                            (pg(1:(end-1))-gamma_mean);
    %prior density for coefficient of basis functions
    prior_tau = reshape([repmat(sigmasqalpha,1,4);...
                        reshape(kron(tau_temp(:,1),ones(nbasis,1)), nbasis, 4)], 4*ngamma,1);
    prior_tau = prior_tau([1:(end-ngamma),(end-(ngamma-2)):end]);  
    log_gamma_prior_curr =  -0.5*(pg(1:(end-1))-zeros(length(prior_tau),1))'*...
                    inv(diag(prior_tau))*(pg(1:(end-1))-zeros(length(prior_tau),1));
    %Loglike at current values
    [loglike_curr]=whittle_like(yobs_tmp,gamma_curr_temp(:,:,1));
    
    log_proposal_curr=log_gamma_curr_temp;
    log_proposal_prop=log_gamma_prop;
    log_target_prop=loglike_prop+log_gamma_prior_prop;
    log_target_curr=loglike_curr+log_gamma_prior_curr;
end

%*************************************************************
%Calculations acceptance probability	
%*************************************************************
PI = min(1,exp(log_target_prop-log_target_curr+log_proposal_curr-log_proposal_prop));
