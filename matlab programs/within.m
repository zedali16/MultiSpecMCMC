function[A,nseg_new,xi_prop,tau_prop,Beta_prop,Phi_prop,seg_temp]=...
        within(chol_index,ts,nexp_temp,tau_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,Phi_temp)
global nobs dimen nBeta M prob_mm1 tmin

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Does the within-model move in the paper
%
%   Input:
%       1) chol_index - s index matrix
%       2) ts - TxN matrix of time series data
%       3) nexp_temp - number of segments
%       5) tau_temp -  vector of smoothing parameters
%       6) xi_curr_temp - current partitions
%       7) nseg_curr_temp - current number of observations in each segment
%       8) Beta_curr_temp - current vector of coefficients
%       9) Phi_temp - which component changed
%       10) Z - together with Phi_temp, indicating how components should be 
%            update
%   Main Outputs:
%       1) PI - acceptance probability
%       2) nseg_new - new number of observations in each segment
%       3) xi_prop - proposed partitions
%       4) tau_prop - proposed smoothing parameters
%       5) Beta_prop - proposed coefficients
%       6) Phi_prop - proposed indicator variable
%
%   Required programs: postBeta1, postBeta2, Beta_derive1, Beta_derive2, whittle_like
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xi_prop = xi_curr_temp;
Beta_prop = Beta_curr_temp;
nseg_new = nseg_curr_temp;
tau_prop = tau_temp;
Phi_prop = Phi_temp;


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
           %log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior)+log(1/2)+log(prob_mm1);
           log_prop_cut_prop=log(1/2)+log(prob_mm1);
        else
           %log_prop_cut_prop=log(1-prob_mm1)-log(nposs_prior)+log(1/3)+log(prob_mm1); 
           log_prop_cut_prop=log(1/3)+log(prob_mm1); 
        end
        if(nseg_new(seg_temp)==tmin || nseg_new(seg_temp+1)==tmin)
           %log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior)+log(1/2)+log(prob_mm1);
           log_prop_cut_curr=log(1/2)+log(prob_mm1);
        else
           %log_prop_cut_curr=log(1-prob_mm1)-log(nposs_prior)+log(1/3)+log(prob_mm1); 
           log_prop_cut_curr=log(1/3)+log(prob_mm1); 
        end
    end

    need = sum(Beta_curr_temp(:,:,seg_temp) - Beta_curr_temp(:,:,seg_temp+1));
    need = need./need;
    need(isnan(need))=0;
    aa = zeros(2^(dimen^2),1);
    for i=1:2^(dimen^2)
          aa(i)=sum(chol_index(i,:)==need);
    end
    Phi_need = find(aa==dimen^2);       
    select = find(chol_index(Phi_need,:)~=0);
    select_inv = find(chol_index(Phi_need,:)==0);
    
    %==========================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the current values
    %==========================================================================
    loglike_curr = 0;
    for j=seg_temp:seg_temp+1 
            if j>1
                yobs_tmp = ts((xi_curr_temp(j-1)+1):xi_curr_temp(j),:);
            else
                yobs_tmp = ts(1:xi_curr_temp(j),:);
            end
         %Loglikelihood at current values
         [log_curr_spec_dens] = whittle_like(yobs_tmp,Beta_curr_temp(:,:,j));  
         loglike_curr = loglike_curr + log_curr_spec_dens;
    end
    
    %=================================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the proposed values
    %=================================================================================
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For coefficient of basis functions that are the same across two segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Phi_need ~=2^(dimen^2)        
        [Beta_out, m_out, m, ~] = Hamilt2(chol_index, Phi_need, seg_temp,...
                                   ts, tau_prop(:,seg_temp), tau_prop(:,seg_temp+1),...
                                   Beta_prop(:,:,seg_temp), Beta_prop(:,:,seg_temp+1), xi_prop, nseg_new(seg_temp));                       
                               
         Beta_prop(:,select_inv,seg_temp) = reshape(Beta_out,nBeta,length(select_inv)); 
         Beta_prop(:,select_inv,seg_temp+1) = Beta_prop(:,select_inv,seg_temp);
         m_curr_1 = -0.5*m'*((1/M)*eye(length(m)))*m;   
         m_prop_1 = -0.5*m_out'*((1/M)*eye(length(m_out)))*m_out;
    else
         m_curr_1 = 0;
         m_prop_1 = 0;
    end    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %For coefficient of basis functions that are different across two segments
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    loglike_prop = 0;
    yobs_tmp_2 = cell(2,1);
    m_prop_2 = 0;
    m_curr_2 = 0;
    for j=seg_temp:seg_temp+1

        [Beta_out, m_out, m, yobs_tmp]=Hamilt1(chol_index, Phi_need,...
                                        j, ts, tau_prop(:,j), Beta_prop(:,:,j), xi_prop);                            
        if j==seg_temp
            yobs_tmp_2{1} = yobs_tmp;
        else
            yobs_tmp_2{2} = yobs_tmp;
        end                    
        
        Beta_prop(:,select,j) = reshape(Beta_out,nBeta,length(select));       
        m_curr_2 = m_curr_2 - 0.5*m'*((1/M)*eye(length(m)))*m;                      
        m_prop_2 = m_prop_2 - 0.5*m_out'*((1/M)*eye(length(m_out)))*m_out; 
    end   
    
    %Loglikelihood at proposed values
    for j=seg_temp:seg_temp+1
        if j==seg_temp
            [log_curr_spec_dens] = whittle_like(yobs_tmp_2{1},Beta_prop(:,:,j)); 
        else
            [log_curr_spec_dens] = whittle_like(yobs_tmp_2{2},Beta_prop(:,:,j));
        end
        loglike_prop = loglike_prop+log_curr_spec_dens;    
    end

    %proposal density
    log_proposal_curr =  log_prop_cut_curr;
    log_proposal_prop =  log_prop_cut_prop;
    
    %target density
    log_prior_cut_prop=0;
	log_prior_cut_curr=0;
    for k=1:nexp_temp-1
        if k==1
            log_prior_cut_prop=-log(nobs-(nexp_temp-k+1)*tmin+1);
			log_prior_cut_curr=-log(nobs-(nexp_temp-k+1)*tmin+1);
		else
			log_prior_cut_prop=log_prior_cut_prop - log(nobs-xi_prop(k-1)-(nexp_temp-k+1)*tmin+1);
			log_prior_cut_curr=log_prior_cut_curr - log(nobs-xi_curr_temp(k-1)-(nexp_temp-k+1)*tmin+1);
        end
    end
    log_target_prop = loglike_prop + log_prior_cut_prop + m_prop_1 + m_prop_2;
    log_target_curr = loglike_curr + log_prior_cut_curr + m_curr_1 + m_curr_2;
    
else    
    
    %*********************************************************
    % If contains only one segment
    %*********************************************************
    nseg_new = nobs;
    seg_temp = 1;
    %=================================================================================
    %Evaluating the Loglikelihood, Priors and Proposals at the proposed values
    %=================================================================================
    Phi_need = 2^(dimen^2);
    
    [Beta_out, m_out, m]=Hamilt1(chol_index, Phi_need,...
                                        1, ts, tau_temp, Beta_curr_temp, xi_curr_temp);
    Beta_prop(:,:,1) = reshape(Beta_out,nBeta,dimen^2);
    m_curr =  - 0.5*m'*((1/M)*eye(length(m)))*m;                      
    m_prop =  - 0.5*m_out'*((1/M)*eye(length(m_out)))*m_out; 

    %Loglike at proposed values
    [loglike_prop] = whittle_like(ts,Beta_prop(:,:,1));
    
    %Loglike at current values
    [loglike_curr] = whittle_like(ts,Beta_curr_temp(:,:,1));
    
    log_target_prop = loglike_prop + m_prop;
    log_target_curr = loglike_curr + m_curr;
    log_proposal_curr =  0;
    log_proposal_prop =  0;
end

%*************************************************************
%Calculations acceptance probability	
%*************************************************************
A = min(1,exp(log_target_prop-log_target_curr+log_proposal_curr-log_proposal_prop));
