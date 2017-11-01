postBeta1 <- function(chol_index,Phi_temp, j, yobs, tau_temp, Beta_temp, xi_temp, global_info){
        
        nbasis <- global_info$nbasis
        nBeta <- global_info$nBeta
        sigmasqalpha <- global_info$sigmasqalpha
        dimen <- global_info$dimen
        if(j>1){
                yobs_tmp <- yobs[(xi_temp[j-1]+1):xi_temp[j],]
        }else{
                yobs_tmp <- yobs[1:xi_temp[j],]
        }
        
        # initial values
         x <- rep(0, chol_index[Phi_temp,]%*%repmat(nBeta,dimen^2,1))
        #x <- as.vector(Beta_temp[,which(chol_index[Phi_temp,]!=0)])

        # optimization
        opt <- trust(Beta_derive1, x, rinit=1, rmax = 10^4, parscale=rep(1,length(x)), 
                     iterlim = 10^4,   fterm = 10^-6, 
                             mterm = 10^-6, minimize = TRUE,  blather = TRUE, 
                     yobs_tmp, chol_index, Phi_temp, tau_temp, Beta_temp, sigmasqalpha, nbasis, global_info)
        Beta_mean <- opt$argument 
        Beta_var <-  matpower(opt$hessian,-1)
        #Beta_var <- solve(opt$hessian)
        list(Beta_mean=Beta_mean, Beta_var=Beta_var, yobs_tmp=yobs_tmp)
}