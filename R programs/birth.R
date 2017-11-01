birth <- function(chol_index,ts,nexp_curr,nexp_prop,
                  tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,Phi_curr_temp,
                  log_move_curr,log_move_prop, global_info){
        
        nobs <- global_info$nobs
        dimen <- global_info$dimen
        nbasis <- global_info$nbasis
        nBeta <- global_info$nBeta
        sigmasqalpha <- global_info$sigmasqalpha
        tmin <- global_info$tmin
        tau_up_limit <- global_info$tau_up_limit
        nBeta <- global_info$nBeta
        
        Beta_prop <- array(0,dim=c(nBeta, dimen^2, nexp_prop))
        tau_prop <- matrix(1, dimen^2, nexp_prop)
        nseg_prop <- matrix(0, nexp_prop, 1)
        xi_prop <- matrix(0, nexp_prop, 1)
        Phi_prop <- matrix(0, nexp_prop,1)
        
        #================================
        # draw segment to spilit
        #================================
        kk <- which(nseg_curr_temp > (2*tmin))
        nposs_seg <- length(kk)
        seg_cut <- kk[sample(1:nposs_seg, 1, replace = TRUE)]
        nposs_cut <- nseg_curr_temp[seg_cut] - 2 * tmin + 1
        #=============================================
        # propose new parameters, Beta, Phi, tau, xi
        #=============================================
        for (jj in 1:nexp_curr){
                if (jj < seg_cut) {
                        xi_prop[jj] <- xi_curr_temp[jj]
                        tau_prop[,jj] <- tau_curr_temp[,jj]
                        nseg_prop[jj] <- nseg_curr_temp[jj]
                        Beta_prop[,,jj] <-  Beta_curr_temp[,,jj]
                        Phi_prop[jj] <- Phi_curr_temp[jj]
                }else if (jj==seg_cut){
                        index <- sample(1:nposs_cut, 1, replace = TRUE)
                        if (seg_cut == 1) {
                                xi_prop[seg_cut] <- index + tmin - 1
                        }else {
                                xi_prop[seg_cut] <- xi_curr_temp[jj-1] - 1 + tmin + index
                        }
                        xi_prop[seg_cut + 1] <- xi_curr_temp[jj]
                        
                        
                        # determin which Cholesky components should change 
                        Phi_prop[seg_cut] <- sample(2:2^(dimen^2), 1, replace = TRUE)
                        Phi_prop[seg_cut + 1] <- Phi_curr_temp[jj]
                        
                        # draw new tausq
                        select <- which(chol_index[Phi_prop[seg_cut],]!=0)
                        select_inv <- which(chol_index[Phi_prop[seg_cut],]==0)
                        zz <- runif(dimen^2) * chol_index[Phi_prop[seg_cut],]
                        uu <- zz/(1-zz)
                        uu[which(uu==0)] <- 1
                        tau_prop[,seg_cut] <- tau_curr_temp[,seg_cut]*uu
                        tau_prop[,seg_cut+1] <- tau_curr_temp[,seg_cut]*(1/uu)
                        
                        # draw bew values for coefficient of basis functions for new birthed segments
                        ff <- 0
                        nseg_prop[seg_cut] <- index + tmin - 1 
                        nseg_prop[seg_cut + 1] <- nseg_curr_temp[jj] - nseg_prop[seg_cut]
                        Phi_need <- Phi_prop[seg_cut]
                        for (k in jj:(jj+1)){
                                if(k==jj){
                                        fit <- postBeta1(chol_index, Phi_need, k, ts, 
                                                         tau_prop[,k], Beta_curr_temp[,,jj], xi_prop, global_info)
                                        Beta_mean_1 <- fit$Beta_mean
                                        Beta_var_1 <- fit$Beta_var
                                        yobs_tmp_1 <- fit$yobs_tmp
                                        if (min(eigen(0.5*(t(fit$Beta_var)+fit$Beta_var))$values)<0){
                                                Beta_prop[,,k] <- Beta_curr_temp[,,jj]
                                                ff <- ff + 1
                                        }else{
                                                Beta_prop[,select,k] <- matrix(rmvnorm(1, fit$Beta_mean, 0.5*(fit$Beta_var + t(fit$Beta_var)),method="chol"),
                                                                               nBeta,length(select))
                                                Beta_prop[,select_inv,k] <- Beta_curr_temp[,select_inv,jj]
                                        }

                                } else{
                                        fit <- postBeta1(chol_index, Phi_need, k, ts, 
                                                         tau_prop[,k], Beta_curr_temp[,,jj], xi_prop, global_info)
                                        Beta_mean_2 <- fit$Beta_mean
                                        Beta_var_2 <- fit$Beta_var
                                        yobs_tmp_2 <- fit$yobs_tmp
                                        if (min(eigen(0.5*(t(fit$Beta_var)+fit$Beta_var))$values)<0){
                                                Beta_prop[,,k] <- Beta_curr_temp[,,jj]
                                                ff <- ff + 1
                                        }else{
                                                Beta_prop[,select,k] <- matrix(rmvnorm(1, fit$Beta_mean, 0.5*(fit$Beta_var + t(fit$Beta_var)),method="chol"),
                                                                               nBeta,length(select))
                                                Beta_prop[,select_inv,k] <- Beta_curr_temp[,select_inv,jj]
                                        }
                                }
                        }
                }else{
                        xi_prop[jj+1] <- xi_curr_temp[jj]
                        tau_prop[,jj+1] <- tau_curr_temp[,jj]
                        nseg_prop[jj+1] <- nseg_curr_temp[jj]
                        Beta_prop[,,jj+1] <- Beta_curr_temp[,,jj]
                        Phi_prop[jj+1] <- Phi_curr_temp[jj]
                }
        }
        
        # calculate Jacobian
        ja <- tau_curr_temp[,seg_cut]/(zz*(1-zz)); ja = ja[ja!=Inf]
        log_jacobian = sum(log(2*ja))
        
        #=========================================================
        # Calculations relatd to proposed values
        #=========================================================
        
        
        #**********************************************************************************
        # Evaluate the Likelihood, proposal and prior densities at the proposed values
        #**********************************************************************************
        log_Beta_prop <- 0
        log_tau_prior_prop <- 0
        log_Beta_prior_prop <- 0
        loglike_prop <- 0
        
        for (jj in seg_cut:(seg_cut+1)){
                if(jj==seg_cut){
                        Beta_mean <- Beta_mean_1
                        Beta_var <- Beta_var_1
                        yobs_tmp <- yobs_tmp_1
                } else{
                        Beta_mean <- Beta_mean_2
                        Beta_var <- Beta_var_2
                        yobs_tmp <- yobs_tmp_2
                }
                pb <- as.vector(Beta_prop[,select,jj]) #, length(Beta_prop[,select,jj]) ,1)
                # proposed density for coefficient of basis function
                log_Beta_prop <- log_Beta_prop - 0.5*t(pb-Beta_mean)%*%matpower(0.5*(Beta_var+t(Beta_var)),-1)%*%(pb-Beta_mean)

                # prior density for coefficient of basis function
                #idx <- cbind( 1*(chol_index[Phi_need, -(1:(dimen + dimen*(dimen-1)/2))] !=1),
                prior_tau <-  matrix(rbind(c(rep(sigmasqalpha,dimen^2-dimen*(dimen-1)/2), tau_prop[-(1:(dimen + dimen*(dimen-1)/2)),jj]),
                                    matrix(kron(tau_prop[,jj],rep(1,nbasis)),nbasis,dimen^2)),nBeta,dimen^2)
                prior_tau <- matrix(prior_tau[,select], length(select)*nBeta,1)
                
                log_Beta_prior_prop <- log_Beta_prior_prop - 0.5*t(pb)%*%matpower(diag(as.vector(prior_tau)),-1)%*%pb

                log_tau_prior_prop <- log_tau_prior_prop - length(select)*log(tau_up_limit)
                log_prop_spec_dens <- whittle_like(yobs_tmp, Beta_prop[,,jj], global_info)
                loglike_prop <- loglike_prop + log_prop_spec_dens
        }
        
        log_seg_prop <- -log(nposs_seg) # Proposal density for segment choice
        log_cut_prop <- -log(nposs_cut) # Proposal density for partition choice
        log_Phi_prop <- -log(2^(dimen^2)-1) # proposal density for component choice
        log_prior_Phi_prop <- -log(2^(dimen^2)) # prior density for component choice
        
        # Evaluating prior density for cut points at proposed values
        log_prior_cut_prop <- 0
        for (k in 1:(nexp_prop-1)){
                if(k==1){
                        log_prior_cut_prop<- - log(nobs-(nexp_prop-k+1)*tmin+1)
                }else{
                        log_prior_cut_prop <- log_prior_cut_prop - log(nobs-xi_prop[k-1]-(nexp_prop-k+1)*tmin+1)
                }
        }
        #Calculating Log Proposal density at Proposed values
        log_proposal_prop <- log_Beta_prop + log_seg_prop + log_move_prop + log_cut_prop + log_Phi_prop
        #Calculating Log Prior density at Proposed values
        log_prior_prop <- log_Beta_prior_prop + log_tau_prior_prop + log_prior_cut_prop + log_prior_Phi_prop
        #Calculating Target density at Proposed values
        log_target_prop <- loglike_prop + log_prior_prop
        
        
        #=========================================================
        # Calculations relatd to current values
        #=========================================================
        
        
        #**********************************************************************************
        # Evaluate the Likelihood, proposal and prior densities at the current values
        #**********************************************************************************
        fit <- postBeta1(chol_index, Phi_need, seg_cut, ts, 
                         tau_curr_temp[,seg_cut], Beta_curr_temp[,,seg_cut], xi_curr_temp, global_info)
        Beta_mean <- fit$Beta_mean
        Beta_var <- fit$Beta_var
        yobs_tmp <- fit$yobs_tmp

        pb <- as.vector(Beta_curr_temp[,select,seg_cut])#, length(Beta_curr_temp[,select,seg_cut]) ,1)
        
        # Current density for coefficient of basis functions
        log_Beta_curr <- -0.5*t(pb-Beta_mean)%*%matpower(0.5*(Beta_var+t(Beta_var)),-1)%*%(pb-Beta_mean)

        # Prior density for coefficient of basis functions at current values 
        prior_tau <-  matrix(rbind(c(rep(sigmasqalpha,dimen^2-dimen*(dimen-1)/2), tau_curr_temp[-(1:(dimen + dimen*(dimen-1)/2)),seg_cut]),
                                   matrix(kron(tau_curr_temp[,seg_cut],rep(1,nbasis)),nbasis,dimen^2)),nBeta,dimen^2)
        prior_tau <- matrix(prior_tau[,select], length(select)*nBeta,1)
        
        log_Beta_prior_curr <-  - 0.5*t(pb)%*%matpower(diag(as.vector(prior_tau)),-1)%*%pb

        
        log_tau_prior_curr <- -length(select)*log(tau_up_limit)
        log_curr_spec_dens <- whittle_like(yobs_tmp, Beta_curr_temp[,,seg_cut], global_info)
        loglike_curr <- log_curr_spec_dens
        
        log_Phi_curr <- -log(1) # proposal for component choice
        #Calculating Log Proposal density at current values
        log_proposal_curr <- log_Beta_curr + log_move_curr + log_Phi_curr
        
        # Evaluating prior density for partition current values
        log_prior_cut_curr <- 0
        for (k in 1:(nexp_curr-1)){
                if (k==1){
                        log_prior_cut_curr <- -log(nobs-(nexp_curr-k+1)*tmin+1)
                }else{
                        log_prior_cut_curr <- log_prior_cut_curr -
                                log(nobs-xi_curr_temp[k-1]-(nexp_curr-k+1)*tmin+1) 
                }
        }
        log_prior_Phi_curr <- -log(2^(dimen^2)) #prior density for component choice
        # Calculating Priors at Current Vlaues
        log_prior_curr <- log_Beta_prior_curr + log_tau_prior_curr + log_prior_cut_curr + log_prior_Phi_curr
        # Evalulating Target densities at current values
        log_target_curr <- loglike_curr + log_prior_curr
        
        #=========================================================
        # Calculations acceptance probability
        #=========================================================
        if(ff>1){
                A <- 0
        }else{
                A <- min(1,exp(log_target_prop - log_target_curr +
                                       log_proposal_curr - log_proposal_prop + log_jacobian))
        }
        list(A=A, nseg_prop=nseg_prop ,xi_prop=xi_prop, 
             tau_prop=tau_prop, Beta_prop=Beta_prop, Phi_prop=Phi_prop)
        
        }