within <- function(chol_index,ts,nexp_temp,tau_temp,
                   xi_curr_temp,nseg_curr_temp,Beta_curr_temp,Phi_temp,global_info){
        nobs <- global_info$nobs
        dimen <- global_info$dimen 
        nBeta <- global_info$nBeta
        M <- global_info$M
        prob_mml <- global_info$prob_mml
        tmin <- global_info$tmin

        
        xi_prop <- xi_curr_temp
        Beta_prop <- Beta_curr_temp
        nseg_new <- nseg_curr_temp
        tau_prop <- tau_temp
        Phi_prop <- Phi_temp
        
        
        if (nexp_temp>1){
                #*********************************************************
                # If contains more than one segments
                #*********************************************************
                
                seg_temp <- sample(1:(nexp_temp-1), 1) # draw segment to cut
                u <- runif(1)
                cut_poss_curr <- xi_curr_temp[seg_temp]
                nposs_prior <- nseg_curr_temp[seg_temp] + nseg_curr_temp[seg_temp+1] - 2*tmin + 1
                
                # Determing if the relocation is big jump or small jump
                if (u < prob_mml) {
                        if (nseg_curr_temp[seg_temp]==tmin&nseg_curr_temp[seg_temp + 1]==tmin){
                                nposs <- 1
                                new_index <- sample(1:nposs, 1, replace = TRUE)
                                cut_poss_new <- xi_curr_temp[seg_temp] - 1 + new_index
                        }else if(nseg_curr_temp[seg_temp]==tmin){
                                nposs <- 2
                                new_index <- sample(1:nposs, 1, replace = TRUE)
                                cut_poss_new <- xi_curr_temp[seg_temp] - 1 + new_index
                        }else if (nseg_curr_temp[seg_temp + 1]==tmin){
                                nposs <- 2
                                new_index <- sample(1:nposs, 1, replace = TRUE)
                                cut_poss_new <- xi_curr_temp[seg_temp] + 1 - new_index
                        }else {
                                nposs <- 3
                                new_index <- sample(1:nposs, 1, replace = TRUE)
                                cut_poss_new <- xi_curr_temp[seg_temp] - 2 + new_index
                        }
                }else{
                        new_index <- sample(1:nposs_prior, 1, replace = TRUE)
                        if (seg_temp>1){
                                cut_poss_new <- sum(nseg_curr_temp[1:(seg_temp - 1)]) - 1 + tmin + new_index
                        }else {
                                cut_poss_new <- -1 + tmin + new_index
                        }
                }
                
                xi_prop[seg_temp] <- cut_poss_new
                if (seg_temp > 1){
                        nseg_new[seg_temp] <- xi_prop[seg_temp] - xi_curr_temp[seg_temp - 1]
                                                                                       
                }else {
                        nseg_new[seg_temp] <- xi_prop[seg_temp]
                }
                nseg_new[seg_temp+1] <- nseg_curr_temp[seg_temp] + nseg_curr_temp[seg_temp+1] - nseg_new[seg_temp]
                
                #=========================================================================================
                # Evaluating the cut Proposal density for the cut-point at the cureent and proposed values
                #=========================================================================================
                if (abs(cut_poss_new - cut_poss_curr) > 1){
                        log_prop_cut_prop <- log(1-prob_mml)-log(nposs_prior)
                        log_prop_cut_curr <- log(1-prob_mml)-log(nposs_prior)
                }else if(nseg_curr_temp[seg_temp]==tmin & nseg_curr_temp[seg_temp+1]==tmin){
                        log_prop_cut_prop <- 0
                        log_prop_cut_curr <- 0
                }else{
                        if(nseg_curr_temp[seg_temp]==tmin||nseg_curr_temp[seg_temp+1]==tmin){
                                log_prop_cut_prop <- log(1/2)+log(prob_mml)
                        }else{
                                log_prop_cut_prop <- log(1/3)+log(prob_mml)
                        }
                        if(nseg_new[seg_temp]==tmin || nseg_new[seg_temp+1]==tmin){
                                log_prop_cut_curr <- log(1/2)+log(prob_mml)
                        }else{
                                log_prop_cut_curr <- log(1/3)+log(prob_mml)
                        }        
                }
                
                need <- colSums(Beta_curr_temp[,,seg_temp] - Beta_curr_temp[,,seg_temp+1])
                need <- need/need
                need[is.na(need)==1] <- 0
                aa <- c()
                for (i in 1:2^(dimen^2)){
                        aa[i] <- sum(chol_index[i,]==need)
                }
                Phi_need <- which(aa==dimen^2)
                select <- which(chol_index[Phi_need,]!=0)
                select_inv <- which(chol_index[Phi_need,]==0)
                
                #==========================================================================
                # Evaluating the Loglikelihood, Priors and Proposals at the current values
                #==========================================================================
                loglike_curr <- 0
                for (j in seg_temp:(seg_temp+1)){
                        if (j>1){
                                yobs_tmp <- ts[(xi_curr_temp[j-1]+1):xi_curr_temp[j],]
                        }else{
                                yobs_tmp <- ts[1:xi_curr_temp[j],]
                        }
                        # loglikelihood at current values
                        log_curr_spec_dens <- whittle_like(yobs_tmp, Beta_curr_temp[,,j], global_info)
                        loglike_curr <- loglike_curr + log_curr_spec_dens
                }
                
                #=================================================================================
                # Evaluating the Loglikelihood, Priors and Proposals at the proposed values
                #=================================================================================
                
                #=================================================================================
                # For coefficient of basis functions that are the same across two segments
                #=================================================================================
                if (Phi_need != 2^(dimen^2)){
                        fit <- Hamilt2(chol_index, Phi_need, seg_temp,
                                       ts, tau_prop[,seg_temp], tau_prop[,seg_temp+1],
                                       Beta_prop[,,seg_temp], Beta_prop[,,seg_temp+1], xi_prop, nseg_new[seg_temp], global_info)
                        
                        Beta_out <- fit$Beta_out; m_out <- fit$m_out; m <- fit$m
                        
                        Beta_prop[,select_inv,seg_temp] <- matrix(Beta_out,nBeta,length(select_inv))
                        Beta_prop[,select_inv,seg_temp+1] <- Beta_prop[,select_inv,seg_temp]
                        m_curr_1 <- -0.5*m%*%((1/M)*eye(length(m)))%*%t(m) 
                        m_prop_1 <- -0.5*m_out%*%((1/M)*eye(length(m_out)))%*%t(m_out)
                }else{
                        m_curr_1 <- 0
                        m_prop_1 <- 0
                }
                
                #=================================================================================
                # For coefficient of basis functions that are different across two segments
                #=================================================================================
                loglike_prop <- 0
                yobs_tmp_2 <- list()
                m_prop_2 <- 0
                m_curr_2 <- 0
                for(j in seg_temp:(seg_temp+1)){
                        
                        fit <- Hamilt1(chol_index, Phi_need, 
                                       j, ts, tau_prop[,j], Beta_prop[,,j], xi_prop, global_info)
                        
                        Beta_out <- fit$Beta_out; m_out <- fit$m_out; m <- fit$m
                        yobs_tmp <- fit$yobs_tmp
                        
                        if (j==seg_temp){
                                yobs_tmp_2[[1]] <- yobs_tmp
                        }else{
                                yobs_tmp_2[[2]] <- yobs_tmp
                        }
                        
                        Beta_prop[,select,j] <- matrix(Beta_out, nBeta, length(select))
                        m_curr_2 = m_curr_2 - 0.5*m%*%((1/M)*eye(length(m)))%*%t(m)                      
                        m_prop_2 = m_prop_2 - 0.5*m_out%*%((1/M)*eye(length(m_out)))%*%t(m_out)
                }
                
                # loglikelihood at proposed values
                for (j in seg_temp:(seg_temp+1)){
                        if (j==seg_temp){
                                log_curr_spec_dens <- whittle_like(yobs_tmp_2[[1]],Beta_prop[,,j],global_info)
                        }else{
                                log_curr_spec_dens <- whittle_like(yobs_tmp_2[[2]],Beta_prop[,,j],global_info)
                        }
                        loglike_prop <- loglike_prop + log_curr_spec_dens 
                }

                
                # proposal density
                log_proposal_curr <-  log_prop_cut_curr
                log_proposal_prop <-  log_prop_cut_prop
                
                # target density
                log_prior_cut_prop <- 0
                log_prior_cut_curr <- 0
                
                for (k in 1:(nexp_temp-1)){
                        if(k==1){
                                log_prior_cut_prop <- -log(nobs-(nexp_temp-k+1)*tmin+1)
                                log_prior_cut_curr <- -log(nobs-(nexp_temp-k+1)*tmin+1)
                        }else{
                                log_prior_cut_prop <- log_prior_cut_prop - log(nobs-xi_prop[k-1]-(nexp_temp-k+1)*tmin+1)
                                log_prior_cut_curr <- log_prior_cut_curr - log(nobs-xi_curr_temp[k-1]-(nexp_temp-k+1)*tmin+1)
                        }
                }
                
                log_target_prop <- loglike_prop + log_prior_cut_prop + m_prop_1 + m_prop_2
                log_target_curr <- loglike_curr + log_prior_cut_curr + m_curr_1 + m_curr_2
        }else{
                #*********************************************************
                # If contains only one segment
                #*********************************************************
                nseg_new <- nobs
                seg_temp <- 1
                #=================================================================================
                # Evaluating the Loglikelihood, Priors and Proposals at the proposed values
                #=================================================================================
                Phi_need <- 2^(dimen^2)
                
                fit <- Hamilt1(chol_index, Phi_need,
                               1, ts, tau_temp, Beta_curr_temp[,,1], xi_curr_temp, global_info)
                Beta_out <- fit$Beta_out; m_out <- fit$m_out; m <- fit$m
                Beta_prop <- matrix(Beta_out, nBeta, dimen^2)
                m_curr <-  - 0.5*m%*%((1/M)*eye(length(m)))%*%t(m)                      
                m_prop <-  - 0.5*m_out%*%((1/M)*eye(length(m_out)))%*%t(m_out) 
                
                # Loglike at proposed values
                loglike_prop <- whittle_like(ts,Beta_prop, global_info)
                
                # loglike at currten values
                loglike_curr <- whittle_like(ts,Beta_curr_temp[,,1], global_info)
                
                log_target_prop <- loglike_prop + m_prop
                log_target_curr <- loglike_curr + m_curr
                log_proposal_curr <-  0
                log_proposal_prop <-  0
        }        
        #*************************************************************
        # Calculations acceptance probability	
        #*************************************************************
                PI <-  min(1,exp(log_target_prop-log_target_curr+
                                         log_proposal_curr-log_proposal_prop))
                
                list(PI=PI,nseg_new=nseg_new,xi_prop=xi_prop,tau_prop=tau_prop,
                     Beta_prop=Beta_prop,Phi_prop=Phi_prop,seg_temp=seg_temp)
}

