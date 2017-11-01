whittle_like <- function(yobs_tmp,Beta,global_info){
        
        dimen <- global_info$dimen
        nBeta <- global_info$nBeta
        Beta_1 <- matrix(0,nBeta,(dimen + dimen*(dimen-1)/2))
        Beta_2 <- matrix(0,nBeta,dimen*(dimen-1)/2)
        Beta_1[,] <- Beta[,1:(dimen + dimen*(dimen-1)/2)]
        Beta_2[,] <- Beta[,-(1:(dimen + dimen*(dimen-1)/2))]
        
        n <- dim(yobs_tmp)[1]
        nfreq <- floor(n/2)
        tt <- (0:nfreq)/(2*nfreq)
        yy <- mvfft(yobs_tmp)/sqrt(n)
        y <- yy[1:(nfreq+1),]
        nf <- dim(y)[1]
        eee <- dim(y)[1]
        
        xx <- lin_basis_func(tt,global_info)
        xx_r <- xx$xx_r; xx_i <- xx$xx_i
        xx_R1 <- xx_r[1,];xx_I1 <- xx_i[1,]
        xx_Re <- xx_r[eee,];xx_Ie <- xx_i[eee,]
        # theta
        theta <- matrix(0, dimen*(dimen-1)/2, nf)
        for(i in 1:(dimen*(dimen-1)/2)){
                theta_real <- xx_r %*% Beta_1[,i+dimen]
                theta_imag <- xx_i %*% Beta_2[,i]
                theta[i,] <- theta_real + 1i*theta_imag 
        }
        
        # delta
        delta_sq <- matrix(0, dimen, nf)
        for(i in 1:dimen){
                delta_sq[i,] <- exp(xx_r%*%Beta_1[,i])
        }
        
        if(dimen==2){
                if(n%%2==1){
                        log_whittle = -sum(log(delta_sq[1,-1]) + log(delta_sq[2,-1]) +
                                        abs(y[-1,1])^2*exp(-xx_r[-1,]%*%Beta_1[,1]) +
                                        abs(y[-1,2] - theta[-1]*y[-1,1])^2*exp(-xx_r[-1,]%*%Beta_1[,2])) -
                                        0.5*(t(log(delta_sq[1,1]) + log(delta_sq[2,1])) + abs(y[1,1])^2*exp(-xx_r[1,]%*%Beta_1[,1]) +
                                        abs(y[1,2] - t(theta[1])*y[1,1])^2*exp(-xx_r[1,]%*%Beta_1[,2]))
                }else{
                        log_whittle = -sum(log(delta_sq[1,2:nfreq]) + log(delta_sq[2,2:nfreq]) +
                                        abs(y[2:nfreq,1])^2*exp(-xx_r[2:nfreq,]%*%Beta_1[,1]) +
                                        abs(y[2:nfreq,2] - theta[2:nfreq]*y[2:nfreq,1])^2*exp(-xx_r[2:nfreq,]%*%Beta_1[,2])) -
                                        0.5*(t(log(delta_sq[1,1]) + log(delta_sq[2,1])) + abs(y[1,1])^2*exp(-xx_r[1,]%*%Beta_1[,1]) +
                                        abs(y[1,2] - t(theta[1])*y[1,1])^2*exp(-xx_r[1,]%*%Beta_1[,2])) -
                                        0.5*(t(log(delta_sq[1,eee]) + log(delta_sq[2,eee])) + abs(y[eee,1])^2*exp(-xx_r[eee,]%*%Beta_1[,1]) +
                                        abs(y[eee,2] - t(theta[eee])*y[eee,1])^2*exp(-xx_r[eee,]%*%Beta_1[,2]))                            
                }
        }else if(dimen==3){
                if(n%%2==1){
                        y1 <- y[2:eee,1]
                        y2 <- y[2:eee,2]
                        y3 <- y[2:eee,3]
                        log_whittle = -sum(log(delta_sq[1,-1]) + log(delta_sq[2,-1]) + log(delta_sq[3,-1]) +
                                        abs(y1)^2*exp(-crossprod(t(xx_r[-1,]),Beta_1[,1])) +
                                        abs(y2 - theta[1,-1]*y1)^2*exp(-crossprod(t(xx_r[-1,]), Beta_1[,2])) +
                                        abs(y3 - theta[2,-1]*y1 - theta[3,-1]*y2)^2*exp(-crossprod(t(xx_r[-1,]), Beta_1[,3])))-
                                        0.5*(t(log(delta_sq[1,1]) + log(delta_sq[2,1]) + log(delta_sq[3,1])) + 
                                        abs(y[1,1])^2*exp(-xx_R1%*%Beta_1[,1]) +
                                        abs(y[1,2] - t(theta[1,1])*y[1,1])^2*exp(-xx_R1%*%Beta_1[,2]) +
                                        abs(y[1,3] - theta[2,1]*y[1,1] - theta[3,1]*y[1,2])^2*exp(-xx_R1%*%Beta_1[,3]))
                }else{
                        y1 <- y[2:nfreq,1]
                        y2 <- y[2:nfreq,2]
                        y3 <- y[2:nfreq,3]
                        log_whittle = -sum(log(delta_sq[1,2:nfreq]) + log(delta_sq[2,2:nfreq]) + log(delta_sq[3,2:nfreq]) +
                                        abs(y1)^2*exp(-crossprod(t(xx_r[2:nfreq,]),Beta_1[,1])) +
                                        abs(y2 - theta[1,2:nfreq]*y1)^2*exp(-crossprod(t(xx_r[2:nfreq,]), Beta_1[,2])) +
                                        abs(y3 - theta[2,2:nfreq]*y1 - theta[3,2:nfreq]*y2)^2*exp(-crossprod(t(xx_r[2:nfreq,]), Beta_1[,3])))-
                                        0.5*(t(log(delta_sq[1,1]) + log(delta_sq[2,1]) + log(delta_sq[3,1])) + 
                                        abs(y[1,1])^2*exp(-xx_R1%*%Beta_1[,1]) +
                                        abs(y[1,2] - t(theta[1,1])*y[1,1])^2*exp(-xx_R1%*%Beta_1[,2]) +
                                        abs(y[1,3] - theta[2,1]*y[1,1] - theta[3,1]*y[1,2])^2*exp(-xx_R1%*%Beta_1[,3])) -
                                        0.5*(t(log(delta_sq[1,eee]) + log(delta_sq[2,eee]) + log(delta_sq[3,eee])) + 
                                        abs(y[eee,1])^2*exp(-xx_Re%*%Beta_1[,1]) +
                                        abs(y[eee,2] - t(theta[1,eee])*y[eee,1])^2*exp(-xx_Re%*%Beta_1[,2]) +
                                        abs(y[eee,3] - theta[2,eee]*y[eee,1] - theta[3,eee]*y[eee,2])^2*exp(-xx_Re%*%Beta_1[,3]))                      
                }
        }
        return(Re(log_whittle))
}


















