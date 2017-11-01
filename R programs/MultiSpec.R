
MultiSpect <- function(zt, nloop, nwarmup, nexp_max, tmin, 
                       prob_mml, nbasis, tau_up_limit, sigmasqalpha,
                       init, nfreq, ee){
        nloop <- nloop
        nwarmup <- nwarmup
        nexp_max <- nexp_max
        
        nbasis <- nbasis
        nBeta <- nbasis + 1
        sigmasqalpha <- sigmasqalpha
        tau_up_limit <- tau_up_limit
        prob_mml <- prob_mml
        tmin <- tmin
        M <- tau_up_limit
        ee <- ee
        
        nfreq_hat <- nfreq
        freq_hat <-(0:nfreq_hat)/(2*nfreq_hat)
        
        
        #====================================
        # Run the estimation procedure
        #====================================
        ts <- zt
        nobs <- dim(ts)[1]
        dimen <- dim(ts)[2]
        global_info <- list(nobs = nobs, nbasis = nbasis, nBeta = nBeta, sigmasqalpha = sigmasqalpha, 
                            prob_mml = prob_mml, tmin = tmin,tau_up_limit = tau_up_limit,dimen=dimen, M=M, ee=ee)
        tausq <- matrix(list(), nexp_max, 1)
        Beta <- matrix(list(), nexp_max, 1)
        spect_hat <- matrix(list(), nexp_max, 1)
        xi <- matrix(list(), nexp_max, 1)
        nseg <- matrix(list(), nexp_max, 1)
        Phi <- matrix(list(), nexp_max, 1)
        tms <- matrix(0, 1, nloop)
        
        
        for (j in 1:nexp_max) {
                tausq[[j]] <- array(1, dim=c(dimen^2, j, nloop + 1))
                Beta[[j]] <- array(0, dim = c(nBeta, dimen^2, j, nloop + 1))
                spect_hat[[j]] <- array(0, dim=c(dimen,dimen,nfreq_hat+1,j,nloop+1))
                xi[[j]] <- matrix(1, j, nloop + 1)
                nseg[[j]] <- matrix(1, j, nloop + 1)
                Phi[[j]] <- matrix(0, j, nloop + 1)
        }
        #=================================
        # initilize the MCMC iteration
        #=================================
        nexp_curr <- rep(nexp_max, nloop+1)
        nexp_curr[1] <- init
        for(j in 1:nexp_curr[1]){
                tausq[[nexp_curr[1]]][, j, 1] <- runif(dimen^2) * tau_up_limit
        }
        
        for(j in 1: nexp_curr[1]){
                if(nexp_curr[1]==1){
                        xi[[nexp_curr[1]]][j, 1] <- nobs
                        nseg[[nexp_curr[1]]][j, 1] <- nobs
                }else{
                        if(j==1){
                                nposs <- nobs-nexp_curr[1]*tmin+1
                                xi[[nexp_curr[1]]][j, 1] <- tmin+sample(1:nposs,1,replace=TRUE)-1
                                nseg[[nexp_curr[1]]][j, 1] <- xi[[nexp_curr[1]]][j,1]
                        } else if(j > 1 & j < nexp_curr[1]){
                                nposs <- nobs-xi[[nexp_curr[1]]][j-1,1]-tmin*(nexp_curr[1]-j+1)+1
                                xi[[nexp_curr[1]]][j,1] <- tmin+sample(1:nposs,1,replace=TRUE)+xi[[nexp_curr[1]]][j-1,1]-1
                                nseg[[nexp_curr[1]]][j,1] <- xi[[nexp_curr[1]]][j,1]-xi[[nexp_curr[1]]][j-1,1]
                        } else{
                                xi[[nexp_curr[1]]][j,1] <- nobs
                                nseg[[nexp_curr[1]]][j,1] <- xi[[nexp_curr[1]]][j, 1] - xi[[nexp_curr[1]]][j-1,1]
                        }
                }
        }
        # index matrix for which compoents of choloskey decomposition changed
        chol_index = chol_ind(dimen)
        
        xi_temp <- xi[[nexp_curr[1]]][, 1]
        tau_temp <- tausq[[nexp_curr[1]]][,,1]
        Beta_temp <- Beta[[nexp_curr[1]]][,,,1]
        Phi[[nexp_curr[1]]][1:(nexp_curr[1]-1),1] <- 2^(dimen^2)
        Phi_temp <- repmat(2^(dimen^2), nexp_curr[1])
        for(j in 1:nexp_curr[1]){
                fit <- postBeta1(chol_index, Phi_temp[j], j, ts, tau_temp[,j], Beta_temp[,,j], xi_temp, global_info)
                Beta[[nexp_curr[1]]][, , j, 1] <- matrix(rmvnorm(1, fit$Beta_mean, fit$Beta_var,method="chol"),nBeta,dimen^2)
        }
        
        # jumping probabilities
        epsilon <- c();met_rat <- c()
        bet_death <- 0;bet_birth <- 0
        with <- 0
        
        #====================================
        # run the loop
        #====================================
        for (p in 1:nloop){
                ptm <- proc.time()[1] 
                if(p%%50==0){
                        cat("iteration",p, "\n")
                }
                #***************************
                # Between Model Move
                #***************************
                kk <- length(which(nseg[[nexp_curr[p]]][,p]>2*tmin))
                
                #*********************************
                # Deciding on birth or death
                #*********************************
                if (kk == 0) {
                        if (nexp_curr[p] == 1) {
                                nexp_prop <- nexp_curr[p]
                                log_move_prop <- 0
                                log_move_curr <- 0
                        }
                        else {
                                nexp_prop <- nexp_curr[p] - 1
                                log_move_prop <- 1
                                if (nexp_prop == 1) {
                                        log_move_curr <- 1
                                }
                                else {
                                        log_move_curr <- log(0.5)
                                }
                        }
                }else {
                        if (nexp_curr[p] == 1) {
                                nexp_prop <- nexp_curr[p] + 1
                                log_move_prop <- 0
                                if (nexp_prop == nexp_max) {
                                        log_move_curr <- 0
                                }
                                else {
                                        log_move_curr <- log(0.5)
                                }
                        }
                        else if (nexp_curr[p] == nexp_max) {
                                nexp_prop <- nexp_curr[p] - 1
                                log_move_prop <- 0
                                if (nexp_prop == 1) {
                                        log_move_curr <- 0
                                }
                                else {
                                        log_move_curr <- log(0.5)
                                }
                        }
                        else {
                                u <- runif(1)
                                if (u < 0.5) {
                                        nexp_prop <- nexp_curr[p] + 1
                                        if (nexp_prop == nexp_max) {
                                                log_move_curr <- 0
                                                log_move_prop <- log(0.5)
                                        }
                                        else {
                                                log_move_curr <- log(0.5)
                                                log_move_prop <- log(0.5)
                                        }
                                }
                                else {
                                        nexp_prop <- nexp_curr[p] - 1
                                        if (nexp_prop == 1) {
                                                log_move_curr <- 0
                                                log_move_prop <- log(0.5)
                                        }
                                        else {
                                                log_move_curr <- log(0.5)
                                                log_move_prop <- log(0.5)
                                        }
                                }
                        }
                }
                xi_curr_temp <- xi[[nexp_curr[p]]][, p] 
                Beta_curr_temp <- array(Beta[[nexp_curr[p]]][,,,p],dim=c(nBeta, dimen^2, nexp_curr[[p]]))
                nseg_curr_temp <- nseg[[nexp_curr[p]]][, p]
                tau_curr_temp <- matrix(tausq[[nexp_curr[p]]][,,p],dimen^2,nexp_curr[p])
                Phi_curr_temp <- Phi[[nexp_curr[p]]][,p]
                
                if(nexp_prop > nexp_curr[p]){
                        # birth
                        birth_move <- birth(chol_index,ts,nexp_curr[p],nexp_prop,
                                            tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,Phi_curr_temp,
                                            log_move_curr,log_move_prop, global_info)
                        met_rat[p] <- birth_move$A
                        nseg_prop <- birth_move$nseg_prop
                        xi_prop <- birth_move$xi_prop
                        tausq_prop <- birth_move$tau_prop
                        Beta_prop <- birth_move$Beta_prop
                        Phi_prop <- birth_move$Phi_prop
                }else if(nexp_prop < nexp_curr[p]){
                        # death
                        death_move <- death(chol_index,ts,nexp_curr[p],nexp_prop,
                                            tau_curr_temp,xi_curr_temp,nseg_curr_temp,Beta_curr_temp,Phi_curr_temp,
                                            log_move_curr,log_move_prop, global_info)
                        met_rat[p] <- death_move$A
                        nseg_prop <- death_move$nseg_prop
                        xi_prop <- death_move$xi_prop
                        tausq_prop <- death_move$tau_prop
                        Beta_prop <- death_move$Beta_prop
                        Phi_prop <- death_move$Phi_prop
                }else{
                        xi_prop <- xi[[nexp_curr[p]]][, p]
                        nseg_prop <- nseg[[nexp_curr[p]]][, p]
                        tausq_prop <- tausq[[nexp_curr[p]]][,,p]
                        Beta_prop <- Beta[[nexp_curr[p]]][,,,p]
                        Phi_prop <- Phi[[nexp_curr[p]]][,p]
                        meta_rat[p] <- 1
                }
                u <- runif(1)
                if (u < met_rat[p]){
                        if (nexp_prop < nexp_curr[p]){
                                bet_death <- bet_death + 1
                        }else if (nexp_prop > nexp_curr[p]){
                                bet_birth <- bet_birth + 1
                        }
                        nexp_curr[p+1] <- nexp_prop
                        xi[[nexp_curr[p+1]]][,p+1] <- xi_prop
                        nseg[[nexp_curr[p+1]]][,p+1] <- nseg_prop
                        tausq[[nexp_curr[p+1]]][,,p+1] <- tausq_prop
                        Beta[[nexp_curr[p+1]]][,,,p+1] <- Beta_prop
                        Phi[[nexp_curr[p+1]]][,p+1] <- Phi_prop
                } else{
                        nexp_curr[p+1] <- nexp_curr[p]
                        xi[[nexp_curr[p+1]]][,p+1] <- xi[[nexp_curr[p+1]]][,p]
                        nseg[[nexp_curr[p+1]]][,p+1] <- nseg[[nexp_curr[p+1]]][,p]
                        tausq[[nexp_curr[p+1]]][,,p+1] <- tausq[[nexp_curr[p+1]]][,,p]
                        Beta[[nexp_curr[p+1]]][,,,p+1] <- Beta[[nexp_curr[p+1]]][,,,p]
                        Phi[[nexp_curr[p+1]]][,p+1] <- Phi[[nexp_curr[p+1]]][,p]
                }
                
                #***************************
                # Within Model Move
                #***************************
                
                # Drawing a new cut point and Beta simultaneously
                # update coeffiecient of linear basis function
                xi_curr_temp <- xi[[nexp_curr[p+1]]][,p+1]
                Beta_curr_temp <- array(Beta[[nexp_curr[p+1]]][,,,p+1],dim=c(nBeta,dimen^2,nexp_curr[[p+1]]))
                tau_curr_temp <- matrix(tausq[[nexp_curr[p+1]]][,,p+1],dimen^2,nexp_curr[p+1])
                nseg_curr_temp <- nseg[[nexp_curr[p+1]]][,p+1]
                Phi_temp <- Phi[[nexp_curr[p+1]]][,p+1] 
                
                within_move <- within(chol_index,ts,nexp_curr[p+1],tau_curr_temp, xi_curr_temp,
                                      nseg_curr_temp, Beta_curr_temp,Phi_temp,global_info)
                epsilon[p] <- within_move$PI
                nseg_new <- within_move$nseg_new
                xi_prop <- within_move$xi_prop
                Beta_prop <- within_move$Beta_prop
                Phi_prop <- within_move$Phi_prop
                seg_temp <- within_move$seg_temp
                
                u <- runif(1)
                if (u < epsilon[p]|| p==1){
                        with <- with + 1
                        if (nexp_curr[p+1]>1){
                                for(j in seg_temp:(seg_temp+1)){
                                        Beta[[nexp_curr[p+1]]][,,j,p+1] <- Beta_prop[,,j]
                                        xi[[nexp_curr[p+1]]][j,p+1] <- xi_prop[j]
                                        nseg[[nexp_curr[p+1]]][j,p+1] <- nseg_new[j]
                                        Phi[[nexp_curr[p+1]]][j,p+1] <- Phi_prop[j]
                                }
                                
                        }else{
                                Beta[[nexp_curr[p+1]]][,,,p+1] <- Beta_prop
                        }
                        
                }else{
                        Beta[[nexp_curr[p+1]]][,,,p+1] <- Beta_curr_temp
                        xi[[nexp_curr[p+1]]][,p+1] <- xi_curr_temp
                        nseg[[nexp_curr[p+1]]][,p+1] <- nseg_curr_temp
                        Phi[[nexp_curr[p+1]]][,p+1] <- Phi_temp
                }
                
                
                # draw tausq
                for (j in 1:nexp_curr[p+1]){
                        for (i in 1:dimen^2){
                                if (is.element(i,1:(dimen + dimen*(dimen-1)/2))){
                                        tau_a <- nbasis/2
                                        tau_b <- sum(Beta[[nexp_curr[p + 1]]][2:nBeta, i, j, p + 1]^2)/2
                                        u <- runif(1)
                                        const1 <- pgamma(1/tau_up_limit, tau_a, tau_b)
                                        const2 <- 1- u * (1-const1)
                                        tausq[[nexp_curr[p+1]]][i,j,p+1] <- 1/qgamma(const2,tau_a,tau_b)
                                }else{
                                        tau_a <- nBeta/2
                                        tau_b <- sum(Beta[[nexp_curr[p + 1]]][1:nBeta, i, j, p + 1]^2)/2
                                        u <- runif(1)
                                        const1 <- pgamma(1/tau_up_limit, tau_a, 1/tau_b)
                                        const2 <-  1- u * (1-const1)
                                        tausq[[nexp_curr[p+1]]][i,j,p+1] <- 1/qgamma(const2,tau_a,tau_b)
                                }
                                
                        }
                }
                
                tms[p] <- proc.time()[1] - ptm 
                cat("run time (sec)",tms[p], "\n")
                #==================================
                # Estimating Spectral Density
                #==================================
                xx <- lin_basis_func(freq_hat,global_info)
                xx_r <- xx$xx_r; xx_i <- xx$xx_i
                
                for (j in 1:nexp_curr[p+1]){
                        g1 <- Beta[[nexp_curr[p+1]]][,1:(dimen + dimen*(dimen-1)/2),j,p+1]
                        g2 <- Beta[[nexp_curr[p+1]]][1:nBeta,-(dimen + dimen*(dimen-1)/2),j,p+1]
                        
                        # theta
                        theta <- matrix(0, dimen*(dimen-1)/2, (nfreq_hat+1))
                        for(i in 1:(dimen*(dimen-1)/2)){
                                theta_real <- xx_r %*% g1[,i+dimen]
                                theta_imag <- xx_i %*% g2[,i]
                                theta[i,] <- theta_real + 1i*theta_imag 
                        }
                        
                        # delta
                        delta_sq <- matrix(0, dimen, (nfreq_hat+1))
                        for(i in 1:dimen){
                                delta_sq[i,] <- exp(xx_r%*%g1[,i])
                        }
                        
                        # produce the spectral density matrix
                        for (k in 1:(nfreq_hat+1)){
                                TT <- diag(dimen)
                                TT[2,1] <- theta[1,k]
                                if (dimen==3){
                                        TT[3,1] <- -theta[2,k]
                                        TT[3,2] <- -theta[3,k]
                                }
                                spect_hat[[nexp_curr[p+1]]][,,k,j,p+1] <-
                                        solve(TT)%*%diag(as.vector(delta_sq[,k]))%*%solve(Conj(t(TT)))
                        }
                        
                }
                
        }
        return(list(spect_hat=spect_hat, freq_hat=freq_hat, tausq=tausq, Beta=Beta, xi=xi, nseg=nseg, 
                    nexp_curr=nexp_curr, Phi=Phi, epsilon=epsilon, bet_birth=bet_birth, bet_death=bet_death,
                    with=with, tms=tms,met_rat=met_rat,nloop=nloop,nwarmup=nwarmup,nfreq_hat=nfreq_hat,
                    nexp_max=nexp_max))
}






