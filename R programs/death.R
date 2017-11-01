death <- function(chol_index,ts,nexp_curr,nexp_prop,
                  tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,Phi_curr_temp,
                  log_move_curr,log_move_prop, global_info){
        
        nobs <- global_info$nobs
        dimen <- global_info$dimen
        nbasis <- global_info$nbasis
        nBeta <- global_info$nBeta
        sigmasqalpha <- global_info$sigmasqalpha
        tmin <- global_info$tmin
        tau_up_limit <- global_info$tau_up_limit
        
        Beta_prop <- array(0,dim=c(nBeta, dimen^2, nexp_prop))
        tau_prop <- matrix(1, dimen^2, nexp_prop)
        nseg_prop <- matrix(0, nexp_prop, 1)
        xi_prop <- matrix(0, nexp_prop, 1)
        Phi_prop <- matrix(0, nexp_prop,1)
        
        #==============================
        # draw a partition to delete
        #==============================
        
        cut_del <- sample(1:(nexp_curr - 1), 1, replace = TRUE)
        
        j <- 0
        for (k in 1:nexp_prop){
                j <- j+1
                if (k==cut_del){
                        #*************************************************************
                        # Calculations related to proposed values	
                        #*************************************************************
                        xi_prop[k] <- xi_curr_temp[j + 1]
                        tau_prop[,k] <- sqrt(tau_curr_temp[,j]*tau_curr_temp[,j+1])
                        nseg_prop[k] <- nseg_curr_temp[j] + nseg_curr_temp[j+1]
                        Phi_prop[k] <- Phi_curr_temp[j+1]
                        
                        #===================================================================
                        # Evaluate the Likelihood at proposed values 
                        #===================================================================
                        
                        need <- colSums(Beta_curr_temp[,,k] - Beta_curr_temp[,,k+1])
                        need <- need/need
                        need[is.na(need)==1] <- 0
                        aa <- c()
                        for (i in 1:2^(dimen^2)){
                                aa[i] <- sum(chol_index[i,]==need)
                        }
                        Phi_need <- which(aa==dimen^2)
                        select <- which(chol_index[Phi_need,]!=0)
                        select_inv <- which(chol_index[Phi_need,]==0)

                        ff <- 0 
                        # Compute mean and variances for coefficents of basis functions
                        fit <- postBeta1(chol_index, Phi_need, k, ts, 
                                         tau_prop[,k], Beta_curr_temp[,,k], xi_prop, global_info)
                        Beta_mean <- fit$Beta_mean
                        Beta_var <- fit$Beta_var
                        yobs_tmp <- fit$yobs_tmp
                        if (min(eigen(0.5*(t(fit$Beta_var)+fit$Beta_var))$values)<0){
                                Beta_prop[,,k] <- Beta_curr_temp[,,k]
                                ff <- ff+1
                        }else{
                                Beta_prop[,select,k] <- matrix(rmvnorm(1, fit$Beta_mean, 0.5*(t(Beta_var)+Beta_var), method="chol"), nBeta,length(select))
                                Beta_prop[,select_inv,k] <- Beta_curr_temp[,select_inv,k]
                        }

                        # loglikelihood at proposed values
                        loglike_prop <- whittle_like(yobs_tmp, Beta_prop[,,k],global_info)
                        
                        
                        #=============================================================================
                        # Evaluate the Proposal Densities at the Proposed values for tau, Phi, and Beta
                        #=============================================================================
                        
                        pb <- as.vector(Beta_prop[,select,k])#, length(Beta_prop[,select,k]) ,1)
                        
                        # Proposed density for coefficient of basis functions
                        log_Beta_prop <-  -0.5*t(pb-Beta_mean)%*%solve(0.5*(Beta_var+t(Beta_var)))%*%(pb-Beta_mean)
                        #log_Beta_prop <- dmvnorm(pb, Beta_mean, Beta_var, log = TRUE)
                        log_seg_prop <- -log(nexp_curr-1) # proposal for segment choice
                        log_Phi_prop <- -log(1) # proposal for component choice
                        
                        # calculate jacobian
                        log_jacobian <- -sum(log(2*(sqrt(tau_curr_temp[select,j]) + sqrt(tau_curr_temp[select,j+1]))^2))
                        
                        # log proposal probability
                        log_proposal_prop <- log_Beta_prop + log_seg_prop + log_move_prop + log_Phi_prop
                        
                        #===========================================================================
                        # Evaluate the PRIOR Densities at the Proposed values for tau, Phi, and Beta
                        #===========================================================================
                        prior_tau <-  matrix(rbind(c(rep(sigmasqalpha,dimen^2-dimen*(dimen-1)/2), tau_prop[-(1:(dimen + dimen*(dimen-1)/2)),k]),
                                                   matrix(kron(tau_prop[,k],rep(1,nbasis)),nbasis,dimen^2)),nBeta,dimen^2)
                        prior_tau <- matrix(prior_tau[,select], length(select)*nBeta,1)
                        
                        log_Beta_prior_prop <- -0.5*t(pb)%*%matpower(diag(as.vector(prior_tau)),-1)%*%pb
  
                        #Prior Density of tausq
                        log_tau_prior_prop <- -length(select)*log(tau_up_limit)
                        #Prior Density of Phi
                        log_Phi_prior_prop <- -log(2^(dimen^2))   
                        log_prior_prop <- log_tau_prior_prop + log_Beta_prior_prop + log_Phi_prior_prop
                        
                        #*************************************************************
                        # Calculations related to current values			
                        #*************************************************************
                                
                        #==================================================================================
                        # Evaluate the Likelihood, Proposal and Prior Densities at the Current values
                        #==================================================================================
                        log_Beta_curr <- 0
                        log_tau_prior_curr <- 0
                        log_Beta_prior_curr <- 0 
                        loglike_curr <- 0
                        for (jj in j:(j+1)){
                                fit <- postBeta1(chol_index, Phi_need, jj, ts, 
                                                 tau_curr_temp[,jj], Beta_curr_temp[,,jj], xi_curr_temp, global_info)
                                Beta_mean <- fit$Beta_mean
                                Beta_var <- fit$Beta_var
                                yobs_tmp <- fit$yobs_tmp
                                
                                pb <- as.vector(Beta_curr_temp[,select,jj])#, length(Beta_curr_temp[,select,jj]) ,1)
                                
                                # Current density for coefficient of basis functions
                                log_Beta_curr <- log_Beta_curr - 
                                                        0.5*t(pb-Beta_mean)%*%solve(0.5*(Beta_var+t(Beta_var)))%*%(pb-Beta_mean)
        
                                prior_tau <-  matrix(rbind(c(rep(sigmasqalpha,dimen^2-dimen*(dimen-1)/2), tau_curr_temp[-(1:(dimen + dimen*(dimen-1)/2)),jj]),
                                                           matrix(kron(tau_curr_temp[,jj],rep(1,nbasis)),nbasis,dimen^2)),nBeta,dimen^2)
                                prior_tau <- matrix(prior_tau[,select], length(select)*nBeta,1)
                                
                                # Prior density for coefficient of basis functions at current values   
                                log_Beta_prior_curr <- log_Beta_prior_curr -
                                                        0.5*(t(pb)%*%matpower(diag(as.vector(prior_tau)),-1)%*%(pb))
      
                                log_curr_spec_dens <- whittle_like(yobs_tmp, Beta_curr_temp[,,jj],global_info)
                                #Loglikelihood at proposed values
                                loglike_curr <- loglike_curr + log_curr_spec_dens
                                
                                # prior density for smoothing parameters
                                log_tau_prior_curr <- log_tau_prior_curr - length(select)*log(tau_up_limit)
                        }
                        
                        log_Phi_curr <- -log(2^(dimen^2)-1) # proposal for component choice
                        log_Phi_prior_curr <- -log(2^(dimen^2)) # prior for component choice
                        
                        # Calculate Log proposal density at current values
                        log_proposal_curr <- log_move_curr + log_Beta_curr + log_Phi_curr
                        
                        # Calculate Priors at Current Vlaues
                        log_prior_curr <- log_Beta_prior_curr + log_tau_prior_curr + log_Phi_prior_curr
                        j <- j+1
                }else{
                        xi_prop[k] <- xi_curr_temp[j]
                        tau_prop[,k] <- tau_curr_temp[,j]
                        nseg_prop[k] <- nseg_curr_temp[j]
                        Beta_prop[,,k] <- Beta_curr_temp[,,j]
                        Phi_prop[k] <- Phi_curr_temp[j]
                }
        }
        
        #=================================================
        # Evaluate Target density at proposed values
        #=================================================
        log_prior_cut_prop <- 0
        for (k in 1:(nexp_prop-1)){
                if (k==1){
                        log_prior_cut_prop <- -log(nobs-(nexp_prop-k+1)*tmin+1)
                }else{
                        log_prior_cut_prop <- log_prior_cut_prop - log(nobs-xi_prop[k-1]-(nexp_prop-k+1)*tmin+1)
                }
        }
        log_target_prop <- loglike_prop + log_prior_prop + log_prior_cut_prop
        
        #==================================================
        # Evaluate Target density at current values
        #==================================================
        log_prior_cut_curr <- 0
        for (k in 1:(nexp_curr-1)){
                if (k==1){
                        log_prior_cut_curr <- -log(nobs-(nexp_curr-k+1)*tmin+1)
                }else{
                        log_prior_cut_curr <- log_prior_cut_curr - log(nobs- xi_curr_temp[k-1]-(nexp_curr-k+1)*tmin+1)       
                }
        }
        log_target_curr <- loglike_curr + log_prior_curr + log_prior_cut_curr
        
        
        #*************************************************************
        # Calculations acceptance probability	
        #*************************************************************
        if(ff>1){
                A <- 0
        }else{
                A = min(1,exp(log_target_prop - log_target_curr + 
                                      log_proposal_curr - log_proposal_prop + log_jacobian))  
        } 
        
        list(A=A, nseg_prop=nseg_prop, xi_prop=xi_prop, 
             tau_prop=tau_prop, Beta_prop=Beta_prop, Phi_prop=Phi_prop)
}