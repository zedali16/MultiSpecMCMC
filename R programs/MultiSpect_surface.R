MultiSpect_surface <- function(zt,output){
        attach(output)
        spectra <- list()
        coh <- list()
        nobs = dim(zt)[1]; dimen=dim(zt)[2]
        if(dimen==2){
                s_11 <- matrix(0,nfreq_hat+1,nobs)
                s_22 <- matrix(0,nfreq_hat+1,nobs)
                coh_21 <- matrix(0,nfreq_hat+1,nobs)
                for (p in 1:nloop){
                        if (p>nwarmup){
                                xi_curr <- xi[[nexp_curr[p]]][,p]
                                spec_hat_curr_11 <- matrix(spect_hat[[nexp_curr[p]]][1,1,,,p],nfreq_hat+1,nexp_curr[[p]])
                                spec_hat_curr_22 <- matrix(spect_hat[[nexp_curr[p]]][2,2,,,p],nfreq_hat+1,nexp_curr[[p]])
                                spec_hat_curr_21 <- abs(matrix(spect_hat[[nexp_curr[p]]][2,1,,,p],nfreq_hat+1,nexp_curr[[p]]))^2/
                                        (matrix(spect_hat[[nexp_curr[p]]][1,1,,,p],nfreq_hat+1,nexp_curr[[p]])*
                                                 matrix(spect_hat[[nexp_curr[p]]][2,2,,,p],nfreq_hat+1,nexp_curr[[p]]))
                                for (j in 1:nexp_curr[p]){
                                        if (j==1){
                                                s_11[,1:xi_curr[j]] <- s_11[,1:xi_curr[j]] + (repmat(spec_hat_curr_11[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)   
                                                s_22[,1:xi_curr[j]] <- s_22[,1:xi_curr[j]] + (repmat(spec_hat_curr_22[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)
                                                coh_21[,1:xi_curr[j]] <- coh_21[,1:xi_curr[j]] + (repmat(spec_hat_curr_21[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)
                                        }else{
                                                s_11[,(xi_curr[j-1]+1):xi_curr[j]]  <- s_11[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_11[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                                s_22[,(xi_curr[j-1]+1):xi_curr[j]]  <- s_22[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_22[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                                coh_21[,(xi_curr[j-1]+1):xi_curr[j]]  <- coh_21[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_21[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                        }
                                }
                        }
                }
                spectra[[1]] <- s_11; spectra[[2]] <- s_22
                coh[[1]] <- coh_21
        }else if (dimen==3){
                s_11 <- matrix(0,nfreq_hat+1,nobs)
                s_22 <- matrix(0,nfreq_hat+1,nobs)
                s_33 <- matrix(0,nfreq_hat+1,nobs)
                coh_21 <- matrix(0,nfreq_hat+1,nobs)
                coh_31 <- matrix(0,nfreq_hat+1,nobs)
                coh_32 <- matrix(0,nfreq_hat+1,nobs)
                for (p in 1:nloop){
                        if (p>nwarmup){
                                xi_curr <- xi[[nexp_curr[p]]][,p]
                                spec_hat_curr_11 <- matrix(spect_hat[[nexp_curr[p]]][1,1,,,p],nfreq_hat+1,nexp_curr[[p]])
                                spec_hat_curr_22 <- matrix(spect_hat[[nexp_curr[p]]][2,2,,,p],nfreq_hat+1,nexp_curr[[p]])
                                spec_hat_curr_33 <- matrix(spect_hat[[nexp_curr[p]]][3,3,,,p],nfreq_hat+1,nexp_curr[[p]])
                                spec_hat_curr_21 <- abs(matrix(spect_hat[[nexp_curr[p]]][2,1,,,p],nfreq_hat+1,nexp_curr[[p]]))^2/
                                        (matrix(spect_hat[[nexp_curr[p]]][1,1,,,p],nfreq_hat+1,nexp_curr[[p]])*
                                                 matrix(spect_hat[[nexp_curr[p]]][2,2,,,p],nfreq_hat+1,nexp_curr[[p]]))
                                spec_hat_curr_31 <- abs(matrix(spect_hat[[nexp_curr[p]]][3,1,,,p],nfreq_hat+1,nexp_curr[[p]]))^2/
                                        (matrix(spect_hat[[nexp_curr[p]]][1,1,,,p],nfreq_hat+1,nexp_curr[[p]])*
                                                 matrix(spect_hat[[nexp_curr[p]]][3,3,,,p],nfreq_hat+1,nexp_curr[[p]]))
                                spec_hat_curr_32 <- abs(matrix(spect_hat[[nexp_curr[p]]][3,2,,,p],nfreq_hat+1,nexp_curr[[p]]))^2/
                                        (matrix(spect_hat[[nexp_curr[p]]][3,3,,,p],nfreq_hat+1,nexp_curr[[p]])*
                                                 matrix(spect_hat[[nexp_curr[p]]][2,2,,,p],nfreq_hat+1,nexp_curr[[p]]))
                                for (j in 1:nexp_curr[p]){
                                        if (j==1){
                                                s_11[,1:xi_curr[j]] <- s_11[,1:xi_curr[j]] + (repmat(spec_hat_curr_11[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)   
                                                s_22[,1:xi_curr[j]] <- s_22[,1:xi_curr[j]] + (repmat(spec_hat_curr_22[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)
                                                s_33[,1:xi_curr[j]] <- s_33[,1:xi_curr[j]] + (repmat(spec_hat_curr_33[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)
                                                coh_21[,1:xi_curr[j]] <- coh_21[,1:xi_curr[j]] + (repmat(spec_hat_curr_21[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)
                                                coh_31[,1:xi_curr[j]] <- coh_31[,1:xi_curr[j]] + (repmat(spec_hat_curr_31[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)
                                                coh_32[,1:xi_curr[j]] <- coh_32[,1:xi_curr[j]] + (repmat(spec_hat_curr_32[,j,drop=FALSE],1,xi_curr[j]))/(nloop - nwarmup)
                                        }else{
                                                s_11[,(xi_curr[j-1]+1):xi_curr[j]]  <- s_11[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_11[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                                s_22[,(xi_curr[j-1]+1):xi_curr[j]]  <- s_22[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_22[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                                s_33[,(xi_curr[j-1]+1):xi_curr[j]]  <- s_33[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_33[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                                coh_21[,(xi_curr[j-1]+1):xi_curr[j]]  <- coh_21[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_21[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                                coh_31[,(xi_curr[j-1]+1):xi_curr[j]]  <- coh_31[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_31[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                                coh_32[,(xi_curr[j-1]+1):xi_curr[j]]  <- coh_32[,(xi_curr[j-1]+1):xi_curr[j]] + (repmat(spec_hat_curr_32[,j,drop=FALSE],1,xi_curr[j]-xi_curr[j-1]))/(nloop-nwarmup)
                                        }
                                        
                                }
                        }
                }
                spectra[[1]] <- s_11; spectra[[2]] <- s_22; spectra[[3]] <- s_33
                coh[[1]] <- coh_21; coh[[2]] <- coh_31; coh[[3]] <- coh_32
        }
        detach(output)
        return(list=c(spectra=spectra,coh=coh))

}

