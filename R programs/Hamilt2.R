Hamilt2 <- function(chol_index, Phi_temp, j, yobs,
                    tau_temp_1, tau_temp_2, Beta_temp_1, Beta_temp_2,
                    xi_temp, nseg, global_info){
        
        nbasis <- global_info$nbasis
        sigmasqalpha <- global_info$sigmasqalpha
        dimen <- global_info$dimen
        nBeta <- global_info$nBeta
        M <- global_info$M
        ee <- global_info$ee
        n <- dim(yobs)[1]
        
        # pick right portion of the data
        if(j==1){
                yobs_tmp <- yobs[1:xi_temp[j+1],]
        }else if(j==(length(xi_temp)-1)){
                yobs_tmp <- yobs[(xi_temp[j-1]+1):n,]
        }else{
                yobs_tmp <- yobs[(xi_temp[j-1]+1):xi_temp[j+1],]
        }
        
        aa <- c()
        for (i in 1:2^(dimen^2)){
                aa[i] <- sum(chol_index[i,]!=chol_index[Phi_temp,])
        }
        Phi_inv <- which(aa==dimen^2)
        select <- which(chol_index[Phi_inv,]!=0)
        
        ll <- nBeta*length(select)
        gr <- Gradient2(yobs_tmp,chol_index, Phi_inv, tau_temp_1, tau_temp_2, 
                        Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg, global_info)
        Beta_old <- matrix(Beta_temp_1[,select],ll,1)
        Beta_out <- Beta_old

        # generate momentum variable
        m <- rmvnorm(1, rep(0, nrow( M*diag(ll))), M*diag(ll))
        # determine leap number and step
        stepsize <- runif(1, 0, 2*ee)
        leap <- sample(1:(2*ceiling(1/ee)), 1, replace = TRUE)
        
        m_out<- m + 0.5*gr*stepsize
        for(i in 1:leap){
                Beta_out <- Beta_out + stepsize*(1/M)*diag(ll)%*%t(m_out)
                Beta_temp_1[,select] <- matrix(Beta_out, nBeta, length(select))
                Beta_temp_2[,select] <- matrix(Beta_out, nBeta, length(select))
                if(i==leap){
                        gr <- Gradient2(yobs_tmp,chol_index, Phi_inv, tau_temp_1, tau_temp_2, 
                                        Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg, global_info)
                        m_out <- m_out + 0.5*gr*stepsize
                }else{
                        gr <- Gradient2(yobs_tmp,chol_index, Phi_inv, tau_temp_1, tau_temp_2, 
                                        Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg, global_info)
                        m_out <- m_out + 1*gr*stepsize
                }
        }
        m_out <- -m_out
        list(Beta_out=Beta_out, m_out=m_out, m=m, yobs_tmp=yobs_tmp)
        
        
}



