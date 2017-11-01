
postBeta2 <- function(chol_index,Phi_temp, j, yobs, 
                      tau_temp_1, tau_temp2, Beta_temp_1, Beta_temp2, 
                      xi_temp, nseg, global_infor){
        
        dimen <- global_info$dimen
        nBeta <- global_info$nBeta
        nbasis <- global_inf$nbasis
        sigmasqalpha <- global_info$sigmasqalpha
        
        # pick right portion of the data
        ee <- dim(yobs)[1]
        if(j==1){
                yobs_tmp <- yobs[1:xi_temp[j+1],] 
        }else if(j==length(xi_temp)-1){
                yobs_tmp <- yobs[(xi_temp[j-1]+1):ee,]
        }else{
                yobs_tmp <- yobs[(xi_temp[j-1]+1):(xi_temp[j-1]+1),]
        }
        
        aa <- matrix(0,2^(dimen^2),1)
        for(i in 1:2^(dimen^2)){
                aa[i] <- sum(chol_index[i,]!=chol_index[Phi_temp,])
        }
        Phi_inv <- which(aa==dimen^2)
        
        # provide initial values
        x <- matrix(0, chol_index[Phi_inv,]%*%repmat(nBeta,dimen^2,1),1)
        
        # optimization
        opt <- trust(Beta_derive2, x, rinit=1, rmax=10^6, parscale=rep(1,length(x)), 
                     iterlim =10^4,  fterm = 10^-6, 
                     mterm = 10^-6, minimize = TRUE,  blather = FALSE,
                     yobs_tmp, chol_index, Phi_inv, tau_temp_1, tau_temp_2, Beta_temp_1, Beta_temp_2, 
                     sigmasqalpha, nbasis, nseg)
        Beta_mean <- opt$argument
        Beta_var <- matpower(opt$hessian,-1)
        
        list(Beta_mean=Beta_mean, Beta_var=Beta_var, yobs_tmp=yobs_tmp)
}