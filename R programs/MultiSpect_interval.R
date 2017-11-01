MultiSpect_interval <- function(zt,output,alphalevel){
        attach(output)
        spectra <- list()
        coh <- list()
        nobs = dim(zt)[1]; dimen=dim(zt)[2]
        if(dimen==2){
                temp1 <- matrix(0,nfreq+1,nloop-nwarmup)
                temp2 <- matrix(0,nfreq+1,nloop-nwarmup)
                temp3 <- matrix(0,nfreq+1,nloop-nwarmup)
                spect11 <- array(0,c(nfreq_hat+1,nobs,3))
                spect22 <- array(0,c(nfreq_hat+1,nobs,3))
                spect22 <- array(0,c(nfreq_hat+1,nobs,3))
                coh21 <- array(0,c(nfreq_hat+1,nobs,3))
                for (j in 1:nobs){
                        if(j%%10==0){
                                cat("completed", j, "of", nobs, "observation", "\n")
                        }
                        for (p in 1:nloop){
                                if(p>nwarmup){
                                        if(length(xi[[nexp_curr[p]]][,p])==1){
                                                k=1
                                        }else{
                                                k = min(which(j <= xi[[nexp_curr[p]]][,p]))
                                        }
                                        temp1[,p-nwarmup] <- spect_hat[[nexp_curr[p]]][1,1,,k,p]
                                        temp2[,p-nwarmup] <- spect_hat[[nexp_curr[p]]][2,2,,k,p]
                                        temp3[,p-nwarmup] <- abs(spect_hat[[nexp_curr[p]]][2,1,,k,p])^2/
                                                (spect_hat[[nexp_curr[p]]][1,1,,k,p]*
                                                         spect_hat[[nexp_curr[p]]][2,2,,k,p])        
                                }
                        }
                        spect11[,j,1] <- apply(temp1,1,mean); spect11[,j,2] = apply(temp1,1,quantile, alphalevel/2); spect11[,j,3] = apply(temp1,1,quantile, 1-alphalevel/2)
                        spect22[,j,1] <- apply(temp2,1,mean); spect22[,j,2] = apply(temp2,1,quantile, alphalevel/2); spect22[,j,3] = apply(temp2,1,quantile, 1-alphalevel/2)
                        coh21[,j,1] <- apply(temp3,1,mean); coh21[,j,2] = apply(temp3,1,quantile, alphalevel/2); coh21[,j,3] = apply(temp3,1,quantile, 1-alphalevel/2)
                }
                return(list(spect11=spect11, spect22=spect22, coh21=coh21))
                
        }else if (dimen==3){
                temp1 <- matrix(0,nfreq+1,nloop-nwarmup)
                temp2 <- matrix(0,nfreq+1,nloop-nwarmup)
                temp3 <- matrix(0,nfreq+1,nloop-nwarmup)
                temp4 <- matrix(0,nfreq+1,nloop-nwarmup)
                temp5 <- matrix(0,nfreq+1,nloop-nwarmup)
                temp6 <- matrix(0,nfreq+1,nloop-nwarmup)
                spect11 <- array(0,c(nfreq_hat+1,nobs,3))
                spect22 <- array(0,c(nfreq_hat+1,nobs,3))
                spect33 <- array(0,c(nfreq_hat+1,nobs,3))
                coh21 <- array(0,c(nfreq_hat+1,nobs,3))
                coh31 <- array(0,c(nfreq_hat+1,nobs,3))
                coh32 <- array(0,c(nfreq_hat+1,nobs,3))
                for (j in 1:nobs){
                        if(j%%10==0){
                                cat("completed", j, "of", nobs, "observation", "\n")
                        }
                        for (p in 1:nloop){
                                if(p>nwarmup){
                                        if(length(xi[[nexp_curr[p]]][,p])==1){
                                                k=1
                                        }else{
                                                k = min(which(j <xi[[nexp_curr[p]]][,p]))
                                        }
                                        temp1[,p-nwarmup] <- spect_hat[[nexp_curr[p]]][1,1,,k,p]
                                        temp2[,p-nwarmup] <- spect_hat[[nexp_curr[p]]][2,2,,k,p]
                                        temp3[,p-nwarmup] <- spect_hat[[nexp_curr[p]]][3,3,,k,p]
                                        temp4[,p-nwarmup] <- abs(spect_hat[[nexp_curr[p]]][2,1,,k,p])^2/
                                                (spect_hat[[nexp_curr[p]]][1,1,,k,p]*
                                                         spect_hat[[nexp_curr[p]]][2,2,,k,p])
                                        temp5[,p-nwarmup] <- abs(spect_hat[[nexp_curr[p]]][3,1,,k,p])^2/
                                                (spect_hat[[nexp_curr[p]]][1,1,,k,p]*
                                                         spect_hat[[nexp_curr[p]]][3,3,,k,p])
                                        temp6[,p-nwarmup] <- abs(spect_hat[[nexp_curr[p]]][3,2,,k,p])^2/
                                                (spect_hat[[nexp_curr[p]]][3,3,,k,p]*
                                                         spect_hat[[nexp_curr[p]]][2,2,,k,p])
                                }

                        }
                        spect11[,j,1] <- apply(temp1,1,mean); spect11[,j,2] = apply(temp1,1,quantile, alphalevel/2); spect11[,j,3] = apply(temp1,1,quantile, 1-alphalevel/2)
                        spect22[,j,1] <- apply(temp2,1,mean); spect22[,j,2] = apply(temp2,1,quantile, alphalevel/2); spect22[,j,3] = apply(temp2,1,quantile, 1-alphalevel/2)
                        spect33[,j,1] <- apply(temp3,1,mean); spect33[,j,2] = apply(temp3,1,quantile, alphalevel/2); spect33[,j,3] = apply(temp3,1,quantile, 1-alphalevel/2)
                        coh21[,j,1] <- apply(temp4,1,mean); coh21[,j,2] = apply(temp4,1,quantile, alphalevel/2); coh21[,j,3] = apply(temp4,1,quantile, 1-alphalevel/2)
                        coh31[,j,1] <- apply(temp5,1,mean); coh31[,j,2] = apply(temp5,1,quantile, alphalevel/2); coh31[,j,3] = apply(temp5,1,quantile, 1-alphalevel/2)
                        coh32[,j,1] <- apply(temp6,1,mean); coh32[,j,2] = apply(temp6,1,quantile, alphalevel/2); coh32[,j,3] = apply(temp6,1,quantile, 1-alphalevel/2)
                }
                return(list(spect11=spect11, spect22=spect22, coh21=coh21))
        }
        detach(output)
}
