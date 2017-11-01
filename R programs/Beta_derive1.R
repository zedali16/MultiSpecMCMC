
Beta_derive1 <- function(x, yobs_tmp, chol_index, Phi_temp, tau_temp,
                         Beta_temp, sigmasqalpha, nbasis, global_info){
        
        dimen <- global_info$dimen
        
        #initilize Beta_1 and Beta_2: Beta_2 is for imaginary components
        nBeta <- nbasis + 1
        Beta_1 <- matrix(0,nBeta,(dimen + dimen*(dimen-1)/2))
        Beta_2 <- matrix(0,nBeta,dimen*(dimen-1)/2)
        
        select <- chol_index[Phi_temp,] * (1:dimen^2)
        select <- select[select!=0]
        
        x <- matrix(x,nBeta,length(select))
        Beta_temp[,select] <- x
        Beta_1[,] <- Beta_temp[,1:(dimen + dimen*(dimen-1)/2)]
        Beta_2[,] <- Beta_temp[,-(1:(dimen + dimen*(dimen-1)/2))]
        
        n <- dim(yobs_tmp)[1]
        nfreq <- floor(n/2); tt <- (0:nfreq)/(2*nfreq)
        yy <- mvfft(yobs_tmp)/sqrt(n); y <- yy[1:(nfreq+1),]; nf <- dim(y)[1]; eee <- dim(y)[1]
        y11 <- y[1,1]; y12 <- y[1,2]; if(dimen==3){y13 <- y[1,3]}
        ye1 <- y[eee,1]; ye2 <- y[eee,2]; if(dimen==3){ye3 <- y[eee,3]}
        xx <- lin_basis_func(tt,global_info); xx_r <- xx$xx_r; xx_i <- xx$xx_i
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
                        y1 <- y[2:eee,1]
                        y2 <- y[2:eee,2]
                        f = -sum(log(delta_sq[1,-1]) + log(delta_sq[2,-1]) +
                                         abs(y[-1,1])^2*exp(-crossprod(t(xx_r[-1,]),Beta_1[,1])) +
                                         abs(y[-1,2] - theta[-1]*y[-1,1])^2*exp(-crossprod(t(xx_r[-1,]), Beta_1[,2]))) -
                                0.5*(t(log(delta_sq[1,1]) + log(delta_sq[2,1])) + abs(y11)^2*exp(-xx_R1%*%Beta_1[,1]) +
                                             abs(y12 - t(theta[1])*y11)^2*exp(-xx_R1%*%Beta_1[,2]))
                }else{
                        y1 <- y[2:nfreq,1]
                        y2 <- y[2:nfreq,2]
                        f = -sum(log(delta_sq[1,2:nfreq]) + log(delta_sq[2,2:nfreq]) +
                                         abs(y1)^2*exp(-crossprod(t(xx_r[2:nfreq,]),Beta_1[,1])) +
                                         abs(y2 - theta[2:nfreq]*y1)^2*exp(-crossprod(t(xx_r[2:nfreq,]),Beta_1[,2]))) -
                                0.5*(t(log(delta_sq[1,1]) + log(delta_sq[2,1])) + abs(y11)^2*exp(-xx_R1%*%Beta_1[,1]) +
                                         abs(y12 - t(theta[1])*y11)^2*exp(-xx_R1%*%Beta_1[,2])) -
                                0.5*(t(log(delta_sq[1,eee]) + log(delta_sq[2,eee])) + abs(ye1)^2*exp(-xx_Re%*%Beta_1[,1]) +
                                         abs(ye2 - t(theta[eee])*ye1)^2*exp(-xx_Re%*%Beta_1[,2]))
                }
                f = f - (0.5*(Beta_1[1,1]%*%t(Beta_1[1,1]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,1])%*%Beta_1[2:nBeta,1])/tau_temp[1])*chol_index[Phi_temp,1] -
                        (0.5*(Beta_1[1,2]%*%t(Beta_1[1,2]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,2])%*%Beta_1[2:nBeta,2])/tau_temp[2])*chol_index[Phi_temp,2] -
                        (0.5*(Beta_1[1,3]%*%t(Beta_1[1,3]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,3])%*%Beta_1[2:nBeta,3])/tau_temp[3])*chol_index[Phi_temp,3] -
                        (0.5*(t(Beta_2[1:nBeta,])%*%Beta_2[1:nBeta,])/tau_temp[4])*chol_index[Phi_temp,4]
                
                gr1 <- matrix(0,nBeta,1);  gr2 <- matrix(0,nBeta,1); gr3 <- matrix(0,nBeta,1); gr4 <- matrix(0,nBeta,1);
                gr1[1] <- Beta_1[1,1]/sigmasqalpha; gr1[2:nBeta] <- Beta_1[2:nBeta,1]/tau_temp[1]
                gr2[1] <- Beta_1[1,2]/sigmasqalpha; gr2[2:nBeta] <- Beta_1[2:nBeta,2]/tau_temp[2]
                gr3[1] <- Beta_1[1,3]/sigmasqalpha; gr3[2:nBeta] <- Beta_1[2:nBeta,3]/tau_temp[3]
                gr4[1:nBeta,1] <- Beta_2[1:nBeta,1]/tau_temp[4]
                h11 <- matrix(0,nBeta,nBeta); h22 <- matrix(0,nBeta,nBeta); h33 <- matrix(0,nBeta,nBeta); h44 <- matrix(0,nBeta,nBeta)
                h42 <- matrix(0,nBeta,nBeta); h32 <- matrix(0,nBeta,nBeta)
                h11[1,1] <- 1/sigmasqalpha; h11[2:nBeta,2:nBeta] <- 1/tau_temp[1]*diag(nbasis)
                h22[1,1] <- 1/sigmasqalpha; h22[2:nBeta,2:nBeta] <- 1/tau_temp[2]*diag(nbasis)
                h33[1,1] <- 1/sigmasqalpha; h33[2:nBeta,2:nBeta] <- 1/tau_temp[3]*diag(nbasis)
                h44[1:nBeta,1:nBeta]  <- 1/tau_temp[4]*diag(nBeta)
                if(n%%2==1){
                        #==============
                        # gradient
                        #==============
                        y1 <- y[2:eee,1];y2 <- y[2:eee,2]
                        xx_R <- xx_r[2:eee,];xx_I <- xx_i[2:eee,]
                        rk <- -y1*Conj(y2) - y2*Conj(y1)
                        ik <- 1i*(-y1*Conj(y2) + y2*Conj(y1))
                        ck <- 2*abs(y1)^2
                        
                        gr1 <- gr1 + crossprod(xx_R, (1-abs(y1)^2*exp(-crossprod(t(xx_R),Beta_1[,1])))) +
                                0.5*(xx_R1*(1-abs(y11)^2*exp(-xx_R1%*%Beta_1[,1])))
                        gr2 <- gr2 + crossprod(xx_R, (1-abs(y2 - theta[2:eee]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*(xx_R1*(1-abs(y12 - t(theta[1])*y11)^2*exp(-xx_R1%*%Beta_1[,2])))
                        temp_mat_31 <- xx_R*rk
                        temp_mat_32 <- ck*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,3])))
                        gr3 <- gr3 + colSums((temp_mat_31 + temp_mat_32)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(-y11*Conj(y12) - y12*Conj(y11))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2)*as.vector(xx_R1%*%Beta_1[,3])*t(xx_R1))
                        temp_mat_41 <- xx_I*ik
                        temp_mat_42 <- ck*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,1])))
                        gr4 <- gr4 + colSums((temp_mat_41 + temp_mat_42)*as.vector(exp(-xx_R%*%Beta_1[,2]))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(1i*(-y11*Conj(y12) + y12*Conj(y11)))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2)*as.vector(xx_I1%*%Beta_2[,1])*t(xx_I1))
                        #==============
                        # hession
                        #==============
                        bigmat_h11 <- kronecker(as.vector(abs(y1)^2*exp(-crossprod(t(xx_R),Beta_1[,1])))*xx_R, t(rep(1,nBeta)))
                        coefmat_h11 <- repmat(xx_R,1,nBeta)
                        h11 <- h11 + matrix(colSums(bigmat_h11*coefmat_h11),nBeta,nBeta) +
                                0.5*(abs(y11)^2*as.vector(exp(-xx_R1%*%Beta_1[,1]))*xx_R1%*%t(xx_R1))
                        
                        bigmat_h22 <- kronecker(as.vector(abs(y2-theta[2:eee]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))
                                                *xx_R,t(rep(1,nBeta)))
                        coefmat_h22 <- repmat(xx_R,1,nBeta)
                        h22 <- h22 + matrix(colSums(bigmat_h22*coefmat_h22),nBeta,nBeta) +
                                0.5*(abs(y12-theta[1]*y11)^2*as.vector(exp(-xx_R1%*%Beta_1[,2]))*xx_R1%*%t(xx_R1))
                        
                        bigmat_h33 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,2]))*ck)*xx_R, t(rep(1,nBeta)))
                        coefmat_h33 <- repmat(xx_R,1,nBeta)
                        h33 <- h33 + matrix(colSums(bigmat_h33*coefmat_h33),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2])*2*abs(y11)^2)*xx_R1%*%t(xx_R1))
                        
                        bigmat_h44 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,2]))*ck)*xx_I, t(rep(1,nBeta)))
                        coefmat_h44 <- repmat(xx_I,1,nBeta)
                        h44 <- h44 + matrix(colSums(bigmat_h44*coefmat_h44),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2])*2*abs(y11)^2)*xx_I1%*%t(xx_I1))
                        
                        bigmat_h42_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck*(crossprod(t(xx_I),Beta_2[,1])))*
                                                          xx_I, t(rep(1,nBeta)))
                        coefmat_h42_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h42_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R), Beta_1[,2]))*ik)*xx_I, t(rep(1,nBeta)))
                        coefmat_h42_2 <- repmat(xx_R, 1,nBeta)
                        h42 <- h42 + t(matrix(colSums(bigmat_h42_1*coefmat_h42_1 
                                                      + bigmat_h42_2*coefmat_h42_2),nBeta,nBeta)) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2*as.vector(xx_I1%*%Beta_2[,1])*xx_I1 +
                                                                                  1i*(-y11*Conj(y12) + y12*Conj(y11))*xx_I1)%*%t(xx_R1))
                        
                        bigmat_h32_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck*(crossprod(t(xx_R),Beta_1[,3])))
                                                  *xx_R, t(rep(1,nBeta)))
                        coefmat_h32_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h32_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*rk)*xx_R, t(rep(1,nBeta)))
                        coefmat_h32_2 <- repmat(xx_R, 1,nBeta)
                        h32 <- h32 + t(matrix(colSums(bigmat_h32_1*coefmat_h32_1 + 
                                                              bigmat_h32_2*coefmat_h32_2),nBeta,nBeta)) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2*as.vector(xx_R1%*%Beta_1[,3])*xx_R1 +
                                                                                  (-y11*Conj(y12) - y12*Conj(y11))*xx_R1)%*%t(xx_R1))
                        h24=t(h42); h23=t(h32)
                }else{
                        #==============
                        # gradient
                        #==============
                        y1 <- y[2:nfreq,1]; y2 <- y[2:nfreq,2]
                        xx_R <- xx_r[2:nfreq,];xx_I <- xx_i[2:nfreq,]
                        rk <- -y1*Conj(y2) - y2*Conj(y1)
                        ik <- 1i*(-y1*Conj(y2) + y2*Conj(y1))
                        ck <- 2*abs(y1)^2
                        
                        gr1 <- gr1 + crossprod(xx_R, (1-abs(y1)^2*exp(-crossprod(t(xx_R),Beta_1[,1])))) +
                                0.5*(xx_R1*(1-abs(y11)^2*exp(-xx_R1%*%Beta_1[,1]))) +
                                0.5*(xx_Re*(1-abs(ye1)^2*exp(-xx_Re%*%Beta_1[,1])))           
                        gr2 <- gr2 + crossprod(xx_R, (1-abs(y2 - theta[2:nfreq]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*(xx_R1*(1-abs(y12 - t(theta[1])*y11)^2*exp(-xx_R1%*%Beta_1[,2]))) +
                                0.5*(xx_Re*(1-abs(ye2 - t(theta[eee])*ye1)^2*exp(-xx_Re%*%Beta_1[,2])))                    
                        temp_mat_31 <- xx_R*rk
                        temp_mat_32 <- ck*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,3])))
                        gr3 <- gr3 + colSums((temp_mat_31 + temp_mat_32)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(-y11*Conj(y12) - y12*Conj(y11))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2)*as.vector(xx_R1%*%Beta_1[,3])*t(xx_R1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(-ye1*Conj(ye2) - ye2*Conj(ye1))*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(ye1)^2)*as.vector(xx_Re%*%Beta_1[,3])*t(xx_Re))
                        temp_mat_41 <- xx_I*ik
                        temp_mat_42 <- ck*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,1])))
                        gr4 <- gr4 + colSums((temp_mat_41 + temp_mat_42)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(1i*(-y11*Conj(y12) + y12*Conj(y11)))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2)*as.vector(xx_I1%*%Beta_2[,1])*t(xx_I1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(1i*(-ye1*Conj(ye2) + ye2*Conj(ye1)))*t(xx_Ie) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(ye1)^2)*as.vector(xx_Ie%*%Beta_2[,1])*t(xx_Ie))
                        #==============
                        # hession
                        #==============
                        bigmat_h11 <- kronecker(as.vector(abs(y1)^2*exp(-crossprod(t(xx_R),Beta_1[,1])))*xx_R, t(rep(1,nBeta)))
                        coefmat_h11 <- repmat(xx_R,1,nBeta)
                        h11 <- h11 + matrix(colSums(bigmat_h11*coefmat_h11),nBeta,nBeta) +
                                0.5*(abs(y11)^2*as.vector(exp(-xx_R1%*%Beta_1[,1]))*xx_R1%*%t(xx_R1)) +
                                0.5*(abs(ye1)^2*as.vector(exp(-xx_Re%*%Beta_1[,1]))*xx_Re%*%t(xx_Re))
                        
                        bigmat_h22 <- kronecker(as.vector(abs(y2-theta[2:nfreq]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))
                                                *xx_R,t(rep(1,nBeta)))
                        coefmat_h22 <- repmat(xx_R,1,nBeta)
                        h22 <- h22 + matrix(colSums(bigmat_h22*coefmat_h22),nBeta,nBeta) +
                                0.5*(abs(y12-theta[1]*y11)^2*as.vector(exp(-xx_R1%*%Beta_1[,2]))*xx_R1%*%t(xx_R1)) +
                                0.5*(abs(ye2-theta[eee]*ye1)^2*as.vector(exp(-xx_Re%*%Beta_1[,2]))*xx_Re%*%t(xx_Re))
                        
                        bigmat_h33 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck)*xx_R, t(rep(1,nBeta)))
                        coefmat_h33 <- repmat(xx_R,1,nBeta)
                        h33 <- h33 + matrix(colSums(bigmat_h33*coefmat_h33),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2])*2*abs(y11)^2)*xx_R1%*%t(xx_R1)) +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,2])*2*abs(ye1)^2)*xx_Re%*%t(xx_Re))
                        
                        bigmat_h44 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck)*xx_I,t(rep(1,nBeta)))
                        coefmat_h44 <- repmat(xx_I,1,nBeta)
                        h44 <- h44 + matrix(colSums(bigmat_h44*coefmat_h44),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2])*2*abs(y11)^2)*xx_I1%*%t(xx_I1)) +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,2])*2*abs(ye1)^2)*xx_Ie%*%t(xx_Ie)) 
                        
                        bigmat_h42_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck*
                                                                    (crossprod(t(xx_I),Beta_2[,1])))*xx_I, t(rep(1,nBeta)))
                        coefmat_h42_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h42_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ik)*xx_I, t(rep(1,nBeta)))
                        coefmat_h42_2 <- repmat(xx_R, 1,nBeta)
                        h42 <- h42 + t(matrix(colSums(bigmat_h42_1*coefmat_h42_1 + 
                                                              bigmat_h42_2*coefmat_h42_2),nBeta,nBeta)) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2*as.vector(xx_I1%*%Beta_2[,1])*xx_I1 +
                                                                                  1i*(-y11*Conj(y12) + y12*Conj(y11))*xx_I1)%*%t(xx_R1)) +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(ye1)^2*as.vector(xx_Ie%*%Beta_2[,1])*xx_Ie +
                                                                                  1i*(-ye1*Conj(ye2) + ye2*Conj(ye1))*xx_Ie)%*%t(xx_Re))
                        
                        bigmat_h32_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck*(crossprod(t(xx_R),Beta_1[,3])))*xx_R, t(rep(1,nBeta)))
                        coefmat_h32_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h32_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*rk)*xx_R, t(rep(1,nBeta)))
                        coefmat_h32_2 <- repmat(xx_R, 1,nBeta)
                        h32 <- h32 + t(matrix(colSums(bigmat_h32_1*coefmat_h32_1 + 
                                                              bigmat_h32_2*coefmat_h32_2),nBeta,nBeta)) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2*as.vector(xx_R1%*%Beta_1[,3])*xx_R1 +
                                                                                  (-y11*Conj(y12) - y12*Conj(y11))*xx_R1)%*%t(xx_R1)) +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(ye1)^2*as.vector(xx_Re%*%Beta_1[,3])*xx_Re +
                                                                                  (-ye1*Conj(ye2) - ye2*Conj(ye1))*xx_Re)%*%t(xx_Re))  
                        h24=t(h42); h23=t(h32)
                }
                h1 <- cbind(h11, matrix(0,nBeta, 2*nBeta+nBeta))
                h2 <- cbind(matrix(0,nBeta,nBeta),h22,-h23,-h24)
                h3 <- cbind(matrix(0,nBeta,nBeta),-h32,h33,matrix(0,nBeta,nBeta))
                h4 <- cbind(matrix(0,nBeta,nBeta),-h42,matrix(0,nBeta,nBeta),h44)
                gr <- c(gr1,gr2,gr3,gr4); h <- rbind(h1,h2,h3,h4); f <- -Re(f)
                gr_index <- (1:(4*nBeta))* kronecker(chol_index[Phi_temp,1:4],rep(1,nBeta))
                gr_index = gr_index[which(gr_index!=0)]
                gr = Re(gr[gr_index]); h = Re(h[gr_index,gr_index])
                
        }else if(dimen==3){
                if(n%%2==1){
                        y1 <- y[2:eee,1]
                        y2 <- y[2:eee,2]
                        y3 <- y[2:eee,3]
                        f = -sum(log(delta_sq[1,-1]) + log(delta_sq[2,-1]) + log(delta_sq[3,-1]) +
                                 abs(y1)^2*exp(-crossprod(t(xx_r[-1,]),Beta_1[,1])) +
                                 abs(y2 - theta[1,-1]*y1)^2*exp(-crossprod(t(xx_r[-1,]), Beta_1[,2])) +
                                 abs(y3 - theta[2,-1]*y1 - theta[3,-1]*y2)^2*exp(-crossprod(t(xx_r[-1,]), Beta_1[,3])))-
                                0.5*(t(log(delta_sq[1,1]) + log(delta_sq[2,1]) + log(delta_sq[3,1])) + 
                                 abs(y11)^2*exp(-xx_R1%*%Beta_1[,1]) +
                                 abs(y12 - t(theta[1,1])*y11)^2*exp(-xx_R1%*%Beta_1[,2]) +
                                 abs(y13 - theta[2,1]*y11 - theta[3,1]*y12)^2*exp(-xx_R1%*%Beta_1[,3]))
                }else{
                        y1 <- y[2:nfreq,1]
                        y2 <- y[2:nfreq,2]
                        y3 <- y[2:nfreq,3]
                        f = -sum(log(delta_sq[1,2:nfreq]) + log(delta_sq[2,2:nfreq]) + log(delta_sq[3,2:nfreq]) +
                                 abs(y1)^2*exp(-crossprod(t(xx_r[2:nfreq,]),Beta_1[,1])) +
                                 abs(y2 - theta[1,2:nfreq]*y1)^2*exp(-crossprod(t(xx_r[2:nfreq,]), Beta_1[,2])) +
                                 abs(y3 - theta[2,2:nfreq]*y1 - theta[3,2:nfreq]*y2)^2*exp(-crossprod(t(xx_r[2:nfreq,]), Beta_1[,3])))-
                                0.5*(t(log(delta_sq[1,1]) + log(delta_sq[2,1]) + log(delta_sq[3,1])) + 
                                 abs(y11)^2*exp(-xx_R1%*%Beta_1[,1]) +
                                 abs(y12 - t(theta[1,1])*y11)^2*exp(-xx_R1%*%Beta_1[,2]) +
                                 abs(y13 - theta[2,1]*y11 - theta[3,1]*y12)^2*exp(-xx_R1%*%Beta_1[,3])) -
                                0.5*(t(log(delta_sq[1,eee]) + log(delta_sq[2,eee]) + log(delta_sq[3,eee])) + 
                                 abs(ye1)^2*exp(-xx_Re%*%Beta_1[,1]) +
                                 abs(ye2 - t(theta[1,eee])*ye1)^2*exp(-xx_Re%*%Beta_1[,2]) +
                                 abs(ye3 - theta[2,eee]*ye1 - theta[3,eee]*ye2)^2*exp(-xx_Re%*%Beta_1[,3]))
                }
                f = f - (0.5*(Beta_1[1,1]%*%t(Beta_1[1,1]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,1])%*%Beta_1[2:nBeta,1])/tau_temp[1])*chol_index[Phi_temp,1] -
                        (0.5*(Beta_1[1,2]%*%t(Beta_1[1,2]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,2])%*%Beta_1[2:nBeta,2])/tau_temp[2])*chol_index[Phi_temp,2] -
                        (0.5*(Beta_1[1,3]%*%t(Beta_1[1,3]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,3])%*%Beta_1[2:nBeta,3])/tau_temp[3])*chol_index[Phi_temp,3] -
                        (0.5*(Beta_1[1,4]%*%t(Beta_1[1,4]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,4])%*%Beta_1[2:nBeta,4])/tau_temp[4])*chol_index[Phi_temp,4] -
                        (0.5*(Beta_1[1,5]%*%t(Beta_1[1,5]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,5])%*%Beta_1[2:nBeta,5])/tau_temp[5])*chol_index[Phi_temp,5] -
                        (0.5*(Beta_1[1,6]%*%t(Beta_1[1,6]))/sigmasqalpha + 0.5*(t(Beta_1[2:nBeta,6])%*%Beta_1[2:nBeta,6])/tau_temp[6])*chol_index[Phi_temp,6] -
                        (0.5*(t(Beta_2[1:nBeta,1])%*%Beta_2[1:nBeta,1])/tau_temp[7])*chol_index[Phi_temp,7] -
                        (0.5*(t(Beta_2[1:nBeta,2])%*%Beta_2[1:nBeta,2])/tau_temp[8])*chol_index[Phi_temp,8] -
                        (0.5*(t(Beta_2[1:nBeta,3])%*%Beta_2[1:nBeta,3])/tau_temp[9])*chol_index[Phi_temp,9] 
        
                gr1 <- matrix(0,nBeta,1);  gr2 <- matrix(0,nBeta,1); gr3 <- matrix(0,nBeta,1); gr4 <- matrix(0,nBeta,1);
                gr5 <- matrix(0,nBeta,1);  gr6 <- matrix(0,nBeta,1); gr7 <- matrix(0,nBeta,1); gr8 <- matrix(0,nBeta,1);
                gr9 <- matrix(0,nBeta,1)                
                gr1[1] <- Beta_1[1,1]/sigmasqalpha; gr1[2:nBeta] <- Beta_1[2:nBeta,1]/tau_temp[1]
                gr2[1] <- Beta_1[1,2]/sigmasqalpha; gr2[2:nBeta] <- Beta_1[2:nBeta,2]/tau_temp[2]
                gr3[1] <- Beta_1[1,3]/sigmasqalpha; gr3[2:nBeta] <- Beta_1[2:nBeta,3]/tau_temp[3]
                gr4[1] <- Beta_1[1,4]/sigmasqalpha; gr1[2:nBeta] <- Beta_1[2:nBeta,4]/tau_temp[4]
                gr5[1] <- Beta_1[1,5]/sigmasqalpha; gr2[2:nBeta] <- Beta_1[2:nBeta,5]/tau_temp[5]
                gr6[1] <- Beta_1[1,6]/sigmasqalpha; gr3[2:nBeta] <- Beta_1[2:nBeta,6]/tau_temp[6]
                gr7[1:nBeta] <- Beta_2[1:nBeta,1]/tau_temp[7]
                gr8[1:nBeta] <- Beta_2[1:nBeta,2]/tau_temp[8]
                gr9[1:nBeta] <- Beta_2[1:nBeta,3]/tau_temp[9]
                
                h11 <- matrix(0,nBeta,nBeta); h22 <- matrix(0,nBeta,nBeta); h33 <- matrix(0,nBeta,nBeta); h44 <- matrix(0,nBeta,nBeta)
                h55 <- matrix(0,nBeta,nBeta); h66 <- matrix(0,nBeta,nBeta); h77 <- matrix(0,nBeta,nBeta); h88 <- matrix(0,nBeta,nBeta)
                h99 <- matrix(0,nBeta,nBeta)
                h42 <- matrix(0,nBeta,nBeta); h53 <- matrix(0,nBeta,nBeta); h56 <- matrix(0,nBeta,nBeta); h59 <- matrix(0,nBeta,nBeta)
                h63 <- matrix(0,nBeta,nBeta); h68 <- matrix(0,nBeta,nBeta); h72 <- matrix(0,nBeta,nBeta); h83 <- matrix(0,nBeta,nBeta)
                h98 <- matrix(0,nBeta,nBeta); h93 <- matrix(0,nBeta,nBeta); 
                h11[1,1] <- 1/sigmasqalpha; h11[2:nBeta,2:nBeta] <- 1/tau_temp[1]*diag(nbasis)
                h22[1,1] <- 1/sigmasqalpha; h22[2:nBeta,2:nBeta] <- 1/tau_temp[2]*diag(nbasis)
                h33[1,1] <- 1/sigmasqalpha; h33[2:nBeta,2:nBeta] <- 1/tau_temp[3]*diag(nbasis)
                h44[1,1] <- 1/sigmasqalpha; h44[2:nBeta,2:nBeta] <- 1/tau_temp[4]*diag(nbasis)
                h55[1,1] <- 1/sigmasqalpha; h55[2:nBeta,2:nBeta] <- 1/tau_temp[5]*diag(nbasis)
                h66[1,1] <- 1/sigmasqalpha; h66[2:nBeta,2:nBeta] <- 1/tau_temp[6]*diag(nbasis)
                h77[1:nBeta,1:nBeta]  <- 1/tau_temp[7]*diag(nBeta)
                h88[1:nBeta,1:nBeta]  <- 1/tau_temp[8]*diag(nBeta)
                h99[1:nBeta,1:nBeta]  <- 1/tau_temp[9]*diag(nBeta)
                
                if(n%%2==1){
                        y1 <- y[2:eee,1];y2 <- y[2:eee,2]; y3 <- y[2:eee,3]
                        xx_R <- xx_r[2:eee,];xx_I <- xx_i[2:eee,]
                        #==============
                        # gradient
                        #==============
                        rk4 <- - y1*Conj(y2) - y2*Conj(y1)
                        ck4 = 2*abs(y1)^2
                        rk5 = -y1*Conj(y3) - y3*Conj(y1)
                        ck5 = 2*abs(y1)^2
                        b = theta[3,]
                        dk5 =  y2*Conj(y1)*b[2:eee] + Conj(y2)*y1*Conj(b[2:eee])
                        rk6 = -y2*Conj(y3) - y3*Conj(y2)
                        ck6 = 2*abs(y2)^2
                        a = theta[2,]
                        dk6 =  y1*Conj(y2)*a[2:eee] + Conj(y1)*y2*Conj(a[2:eee])
                        ik7 = (1i)*(-y1*Conj(y2) + y2*Conj(y1))
                        ck7 = 2*abs(y1)^2
                        ik8 = (1i)*(-y1*Conj(y3) + y3*Conj(y1))
                        ck8 = 2*abs(y1)^2
                        dk8 = (1i)*(-y2*Conj(y1)*b[2:eee] + Conj(y2*Conj(y1)*b[2:eee]))
                        ik9 = (1i)*(-y2*Conj(y3) + y3*Conj(y2));
                        ck9 = 2*abs(y2)^2
                        dk9 = (1i)*(-y1*Conj(y2)*a[2:eee] + Conj(y1*Conj(y2)*a[2:eee]))
                        
                        gr1 <- gr1 + t(xx_R)%*%(1-abs(y1)^2*exp(-xx_R%*%Beta_1[,1])) +
                                0.5*(xx_R1*(1-abs(y11)^2*exp(-xx_R1%*%Beta_1[,1])))
                        gr2 <- gr2 + crossprod(xx_R, (1-abs(y2 - theta[1,2:eee]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*(xx_R1*(1-abs(y12 - t(theta[1,1])*y11)^2*exp(-xx_R1%*%Beta_1[,2])))
                        gr3 <- gr3 + crossprod(xx_R, (1-abs(y3 - theta[2,2:eee]*y1 - theta[3,2:eee]*y2 )^2*exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*(xx_R1*(1-abs(y13 - t(theta[2,1])*y11 - t(theta[3,1])*y12)^2*exp(-xx_R1%*%Beta_1[,3])))
                        temp_mat_41 <- xx_R*rk4
                        temp_mat_42 <- ck4*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,4])))
                        gr4 <- gr4 + colSums((temp_mat_41 + temp_mat_42)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(-y11*Conj(y12) - y12*Conj(y11))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2)*as.vector(xx_R1%*%Beta_1[,4])*t(xx_R1))
                        temp_mat_51 <- xx_R*rk5
                        temp_mat_52 <- ck5*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,5])))
                        temp_mat_53 <- xx_R*dk5
                        gr5 <- gr5 + colSums((temp_mat_51 + temp_mat_52 + temp_mat_53)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y11*Conj(y13) - y13*Conj(y11))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2)*as.vector(xx_R1%*%Beta_1[,5])*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(y12*Conj(y11)*(b[1]) + Conj(y12*Conj(y11)*(b[1])))*t(xx_R1))
                        temp_mat_61 <- xx_R*rk6
                        temp_mat_62 <- ck6*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,6])))
                        temp_mat_63 <- xx_R*dk6
                        gr6 <- gr6 + colSums((temp_mat_61 + temp_mat_62 + temp_mat_63)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y12*Conj(y13) - y13*Conj(y12))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y12)^2)*as.vector(xx_R1%*%Beta_1[,6])*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(y11*Conj(y12)*(a[1]) + Conj(y11*Conj(y12)*(a[1])))*t(xx_R1))
                        temp_mat_71 <- xx_I*ik7
                        temp_mat_72 <- ck7*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,1])))
                        gr7 <- gr7 + colSums((temp_mat_71 + temp_mat_72)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(Im(-y11*Conj(y12) + y12*Conj(y11)))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2)*as.vector(xx_I1%*%Beta_2[,1])*t(xx_I1))
                        temp_mat_81 <- xx_I*ik8
                        temp_mat_82 <- ck8*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,2])))
                        temp_mat_83 <- xx_I*dk8
                        gr8 <- gr8 + colSums((temp_mat_81 + temp_mat_82 + temp_mat_83)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Im(-y11*Conj(y13) + y13*Conj(y11)))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2)*as.vector(xx_I1%*%Beta_2[,2])*t(xx_I1) +
                                              (1i)*as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y11*Conj(y12)*(b[1]) + Conj(y11*Conj(y12)*(b[1])))*t(xx_I1))                        
                        temp_mat_91 <- xx_I*ik9
                        temp_mat_92 <- ck9*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,3])))
                        temp_mat_93 <- xx_I*dk9
                        gr9 <- gr9 + colSums((temp_mat_91 + temp_mat_92 + temp_mat_93)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Im(-y12*Conj(y13) + y13*Conj(y12)))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y12)^2)*as.vector(xx_I1%*%Beta_2[,3])*t(xx_I1) +
                                              (1i)*as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y11*Conj(y12)*(a[1]) + Conj(y11*Conj(y12)*(a[1])))*t(xx_I1)) 
                        #==============
                        # Hession
                        #==============
                        bigmat_h11 <- kronecker(as.vector(abs(y1)^2*exp(-crossprod(t(xx_R),Beta_1[,1])))*xx_R, t(rep(1,nBeta)))
                        coefmat_h11 <- repmat(xx_R,1,nBeta)
                        h11 <- h11 + matrix(colSums(bigmat_h11*coefmat_h11),nBeta,nBeta) +
                                0.5*(abs(y11)^2*as.vector(exp(-xx_R1%*%Beta_1[,1]))*xx_R1%*%t(xx_R1))
                
                        bigmat_h22 <- kronecker(as.vector(abs(y2-theta[1,2:eee]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))
                                                *xx_R,t(rep(1,nBeta)))
                        coefmat_h22 <- repmat(xx_R,1,nBeta)
                        h22 <- h22 + matrix(colSums(bigmat_h22*coefmat_h22),nBeta,nBeta) +
                                0.5*(abs(y12-theta[1,1]*y11)^2*as.vector(exp(-xx_R1%*%Beta_1[,2]))*xx_R1%*%t(xx_R1))
                        
                        bigmat_h33 <- kronecker(as.vector(abs(y3-theta[2,2:eee]*y1 - theta[3,2:eee]*y2)^2*exp(-crossprod(t(xx_R),Beta_1[,3])))
                                                *xx_R,t(rep(1,nBeta)))
                        coefmat_h33 <- repmat(xx_R,1,nBeta)
                        h33 <- h33 + matrix(colSums(bigmat_h33*coefmat_h33),nBeta,nBeta) +       
                                0.5*(abs(y13-theta[2,1]*y11 - theta[3,1]*y12)^2*as.vector(exp(-xx_R1%*%Beta_1[,3]))*xx_R1%*%t(xx_R1))
                        
                        bigmat_h44 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,2]))*ck4)*xx_R, t(rep(1,nBeta)))
                        coefmat_h44 <- repmat(xx_R,1,nBeta)
                        h44 <- h44 + matrix(colSums(bigmat_h44*coefmat_h44),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2])*2*abs(y11)^2)*xx_R1%*%t(xx_R1))
                        
                        bigmat_h55 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,3]))*ck5)*xx_R, t(rep(1,nBeta)))
                        coefmat_h55 <- repmat(xx_R,1,nBeta)
                        h55 <- h55 + matrix(colSums(bigmat_h55*coefmat_h55),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3])*2*abs(y11)^2)*xx_R1%*%t(xx_R1))                        
                        
                        bigmat_h66 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,3]))*ck6)*xx_R, t(rep(1,nBeta)))
                        coefmat_h66 <- repmat(xx_R,1,nBeta)
                        h66 <- h66 + matrix(colSums(bigmat_h66*coefmat_h66),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3])*2*abs(y12)^2)*xx_R1%*%t(xx_R1))                              
                        
                        bigmat_h77 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,2]))*ck7)*xx_I, t(rep(1,nBeta)))
                        coefmat_h77 <- repmat(xx_I,1,nBeta)
                        h77 <- h77 + matrix(colSums(bigmat_h77*coefmat_h77),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2])*2*abs(y11)^2)*xx_I1%*%t(xx_I1))                          
                        
                        bigmat_h88 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,3]))*ck9)*xx_I, t(rep(1,nBeta)))
                        coefmat_h88 <- repmat(xx_I,1,nBeta)
                        h88 <- h88 + matrix(colSums(bigmat_h88*coefmat_h88),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3])*2*abs(y11)^2)*xx_I1%*%t(xx_I1))                          

                        bigmat_h99 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,3]))*ck8)*xx_I, t(rep(1,nBeta)))
                        coefmat_h99 <- repmat(xx_I,1,nBeta)
                        h99 <- h99 + matrix(colSums(bigmat_h99*coefmat_h99),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3])*2*abs(y12)^2)*xx_I1%*%t(xx_I1))                            
                        
                        bigmat_h42_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck4*(crossprod(t(xx_R),Beta_1[,4])))
                                                  *xx_R, t(rep(1,nBeta)))
                        coefmat_h42_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h42_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*rk4)*xx_R, t(rep(1,nBeta)))
                        coefmat_h42_2 <- repmat(xx_R, 1,nBeta)
                        h42 <- h42 + t(matrix(colSums(bigmat_h42_1*coefmat_h42_1 + 
                                        bigmat_h42_2*coefmat_h42_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2*as.vector(xx_R1%*%Beta_1[,4])*xx_R1 +
                                        (-y11*Conj(y12) - y12*Conj(y11))*xx_R1)%*%t(xx_R1))
                       
                        bigmat_h53_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ck5*(crossprod(t(xx_R),Beta_1[,5])))
                                                  *xx_R, t(rep(1,nBeta)))
                        coefmat_h53_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h53_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*rk5)*xx_R, t(rep(1,nBeta)))
                        coefmat_h53_2 <- repmat(xx_R, 1,nBeta)
                        bigmat_h53_3 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*dk5)*xx_R, t(rep(1,nBeta)))
                        coefmat_h53_3 <- repmat(xx_R, 1,nBeta)                        
                        h53 <- h53 + t(matrix(colSums(bigmat_h53_1*coefmat_h53_1 + 
                                        bigmat_h53_2*coefmat_h53_2 +
                                        bigmat_h53_3*coefmat_h53_3),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2*as.vector(xx_R1%*%Beta_1[,5])*xx_R1 +
                                        (-y11*Conj(y13) - y13*Conj(y11))*xx_R1 + 
                                        (y12*Conj(y11)*b[1] + Conj(y12*Conj(y11)*b[1]))*xx_R1)%*%t(xx_R1))
                        
                        bigmat_h56_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y1)*y2)*xx_R, t(rep(1,nBeta)))
                        coefmat_h56_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h56_2 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y2)*y1)*xx_R, t(rep(1,nBeta)))
                        coefmat_h56_2 <- repmat(xx_R, 1,nBeta)
                        h56 <- Re(h56 + t(matrix(colSums(bigmat_h56_1*coefmat_h56_1 + 
                                        bigmat_h56_2*coefmat_h56_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Conj(y11)*y12*xx_R1 +
                                        Conj(y12)*y11*xx_R1)%*%t(xx_R1))) 
                        
                        bigmat_h59_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y1)*y2)*xx_R, t(rep(1,nBeta)))
                        coefmat_h59_1 <- repmat(xx_I,1,nBeta)
                        bigmat_h59_2 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y2)*y1)*xx_R, t(rep(1,nBeta)))
                        coefmat_h59_2 <- repmat(xx_I, 1,nBeta)
                        h59 <- Im(h59 - t(matrix(colSums(bigmat_h59_1*coefmat_h59_1 - 
                                        bigmat_h59_2*coefmat_h59_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Conj(y11)*y12*xx_R1 +
                                        Conj(y12)*y11*xx_R1)%*%t(xx_I1)))  
                        
                        bigmat_h63_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ck6*(crossprod(t(xx_R),Beta_1[,6])))
                                                  *xx_R, t(rep(1,nBeta)))
                        coefmat_h63_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h63_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*rk6)*xx_R, t(rep(1,nBeta)))
                        coefmat_h63_2 <- repmat(xx_R, 1,nBeta)
                        bigmat_h63_3 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*dk6)*xx_R, t(rep(1,nBeta)))
                        coefmat_h63_3 <- repmat(xx_R, 1,nBeta)                        
                        h63 <- h63 + t(matrix(colSums(bigmat_h63_1*coefmat_h63_1 + 
                                        bigmat_h63_2*coefmat_h63_2 +
                                        bigmat_h63_3*coefmat_h63_3),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2*as.vector(xx_R1%*%Beta_1[,6])*xx_R1 +
                                        (-y12*Conj(y13) - y13*Conj(y12))*xx_R1 + 
                                        (y12*Conj(y11)*a[1] + Conj(y12*Conj(y11)*a[1]))*xx_R1)%*%t(xx_R1))  
                        
                        bigmat_h68_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y2)*y1)*xx_R, t(rep(1,nBeta)))
                        coefmat_h68_1 <- repmat(xx_I,1,nBeta)
                        bigmat_h68_2 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y1)*y2)*xx_R, t(rep(1,nBeta)))
                        coefmat_h68_2 <- repmat(xx_I, 1,nBeta)
                        h68 <- Im(h68 - t(matrix(colSums(bigmat_h68_1*coefmat_h68_1 - 
                                        bigmat_h68_2*coefmat_h68_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Conj(y12)*y11*xx_R1 +
                                        Conj(y11)*y12*xx_R1)%*%t(xx_I1)))       
                        
                        bigmat_h72_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck7*(crossprod(t(xx_I),Beta_2[,1])))
                                                  *xx_I, t(rep(1,nBeta)))
                        coefmat_h72_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h72_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ik7)*xx_I, t(rep(1,nBeta)))
                        coefmat_h72_2 <- repmat(xx_R, 1,nBeta)
                        h72 <- h72 + t(matrix(colSums(bigmat_h72_1*coefmat_h72_1 + 
                                        bigmat_h72_2*coefmat_h72_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2*as.vector(xx_I1%*%Beta_2[,1])*xx_I1 +
                                        Im(-y11*Conj(y12) + y12*Conj(y11))*xx_I1)%*%t(xx_R1))

                        bigmat_h83_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ck8*(crossprod(t(xx_I),Beta_2[,2])))
                                                  *xx_I, t(rep(1,nBeta)))
                        coefmat_h83_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h83_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ik8)*xx_I, t(rep(1,nBeta)))
                        coefmat_h83_2 <- repmat(xx_R, 1,nBeta)
                        bigmat_h83_3 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*dk8)*xx_I, t(rep(1,nBeta)))
                        coefmat_h83_3 <- repmat(xx_R, 1,nBeta)                        
                        h83 <- h83 + t(matrix(colSums(bigmat_h83_1*coefmat_h83_1 + 
                                        bigmat_h83_2*coefmat_h83_2 +
                                        bigmat_h83_3*coefmat_h83_3),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2*as.vector(xx_I1%*%Beta_2[,2])*xx_I1 +
                                        (-y11*Conj(y13) + y13*Conj(y11))*xx_I1 + 
                                        Im(-y12*Conj(y11)*b[1] + Conj(y12*Conj(y11)*b[1]))*xx_I1)%*%t(xx_R1))  
                        
                        bigmat_h93_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ck9*(crossprod(t(xx_I),Beta_2[,3])))
                                                  *xx_I, t(rep(1,nBeta)))
                        coefmat_h93_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h93_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ik9)*xx_I, t(rep(1,nBeta)))
                        coefmat_h93_2 <- repmat(xx_R, 1,nBeta)
                        bigmat_h93_3 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*dk9)*xx_I, t(rep(1,nBeta)))
                        coefmat_h93_3 <- repmat(xx_R, 1,nBeta)                        
                        h93 <- h93 + t(matrix(colSums(bigmat_h93_1*coefmat_h93_1 + 
                                        bigmat_h93_2*coefmat_h93_2 +
                                        bigmat_h93_3*coefmat_h93_3),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y12)^2*as.vector(xx_I1%*%Beta_2[,3])*xx_I1 +
                                        (-y12*Conj(y13) + y13*Conj(y12))*xx_I1 + 
                                        Im(-y11*Conj(y12)*a[1] + Conj(y11*Conj(y12)*a[1]))*xx_I1)%*%t(xx_R1))    
                        
                        bigmat_h98_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y1)*y2)*xx_I, t(rep(1,nBeta)))
                        coefmat_h98_1 <- repmat(xx_I,1,nBeta)
                        bigmat_h98_2 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y2)*y1)*xx_I, t(rep(1,nBeta)))
                        coefmat_h98_2 <- repmat(xx_I, 1,nBeta)
                        h98 <- Re(h98 + t(matrix(colSums(bigmat_h98_1*coefmat_h98_1 + 
                                        bigmat_h98_2*coefmat_h98_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Conj(y11)*y12*xx_I1 +
                                        Conj(y12)*y11*xx_I1)%*%t(xx_I1)))  
                        
                        h24=t(h42); h35=t(h53); h65=t(h56); h95=t(h59); h36=t(h63);
                        h86=t(h68); h27=t(h72); h38=t(h83); h39=t(h93); h89=t(h98);
                }else{
                        y1 <- y[2:nfreq,1];y2 <- y[2:nfreq,2];y3 <- y[2:nfreq,3]
                        xx_R <- xx_r[2:nfreq,];xx_I <- xx_i[2:nfreq,]
                        #==============
                        # gradient
                        #==============
                        rk4 <- - y1*Conj(y2) - y2*Conj(y1)
                        ck4 = 2*abs(y1)^2
                        rk5 = -y1*Conj(y3) - y3*Conj(y1)
                        ck5 = 2*abs(y1)^2
                        b = theta[3,]
                        dk5 =  y2*Conj(y1)*b[2:nfreq] + Conj(y2)*y1*Conj(b[2:nfreq])
                        rk6 = -y2*Conj(y3) - y3*Conj(y2)
                        ck6 = 2*abs(y2)^2
                        a = theta[2,]
                        dk6 =  y1*Conj(y2)*a[2:nfreq] + Conj(y1)*y2*Conj(a[2:nfreq])
                        ik7 = (1i)*(-y1*Conj(y2) + y2*Conj(y1))
                        ck7 = 2*abs(y1)^2
                        ik8 = (1i)*(-y1*Conj(y3) + y3*Conj(y1))
                        ck8 = 2*abs(y1)^2
                        dk8 = (1i)*(-y2*Conj(y1)*b[2:nfreq] + Conj(y2*Conj(y1)*b[2:nfreq]))
                        ik9 = (1i)*(-y2*Conj(y3) + y3*Conj(y2));
                        ck9 = 2*abs(y2)^2
                        dk9 = (1i)*(-y1*Conj(y2)*a[2:nfreq] + Conj(y1*Conj(y2)*a[2:nfreq]))
                        
                        
                        gr1 <- gr1 + t(xx_R)%*%(1-abs(y1)^2*exp(-xx_R%*%Beta_1[,1])) +
                                0.5*(xx_R1*(1-abs(y11)^2*exp(-xx_R1%*%Beta_1[,1])))+
                                0.5*(xx_Re*(1-abs(ye1)^2*exp(-xx_Re%*%Beta_1[,1])))
                        gr2 <- gr2 + crossprod(xx_R, (1-abs(y2 - theta[1,2:nfreq]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*(xx_R1*(1-abs(y12 - theta[1,1]*y11)^2*exp(-xx_R1%*%Beta_1[,2])))+
                                0.5*(xx_Re*(1-abs(ye2 - theta[1,eee]*ye1)^2*exp(-xx_Re%*%Beta_1[,2])))
                        gr3 <- gr3 + crossprod(xx_R, (1-abs(y3 - theta[2,2:nfreq]*y1 - theta[3,2:nfreq]*y2 )^2*exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*(xx_R1*(1-abs(y13 - t(theta[2,1])*y11 - t(theta[3,1])*y12)^2*exp(-xx_R1%*%Beta_1[,3]))) +
                                0.5*(xx_Re*(1-abs(ye3 - t(theta[2,eee])*ye1 - t(theta[3,eee])*ye2)^2*exp(-xx_Re%*%Beta_1[,3])))
                        temp_mat_41 <- xx_R*rk4
                        temp_mat_42 <- ck4*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,4])))
                        gr4 <- gr4 + colSums((temp_mat_41 + temp_mat_42)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(-y11*Conj(y12) - y12*Conj(y11))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2)*as.vector(xx_R1%*%Beta_1[,4])*t(xx_R1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(-ye1*Conj(ye2) - ye2*Conj(ye1))*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(ye1)^2)*as.vector(xx_Re%*%Beta_1[,4])*t(xx_Re))
                        temp_mat_51 <- xx_R*rk5
                        temp_mat_52 <- ck5*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,5])))
                        temp_mat_53 <- xx_R*dk5
                        gr5 <- gr5 + colSums((temp_mat_51 + temp_mat_52 + temp_mat_53)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y11*Conj(y13) - y13*Conj(y11))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2)*as.vector(xx_R1%*%Beta_1[,5])*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(y12*Conj(y11)*(b[1]) + Conj(y12*Conj(y11)*(b[1])))*t(xx_R1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(-ye1*Conj(ye3) - ye3*Conj(ye1))*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(ye1)^2)*as.vector(xx_Re%*%Beta_1[,5])*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(ye2*Conj(ye1)*(b[eee]) + Conj(ye2*Conj(ye1)*(b[eee])))*t(xx_Re))
                        temp_mat_61 <- xx_R*rk6
                        temp_mat_62 <- ck6*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,6])))
                        temp_mat_63 <- xx_R*dk6
                        gr6 <- gr6 + colSums((temp_mat_61 + temp_mat_62 + temp_mat_63)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y12*Conj(y13) - y13*Conj(y12))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y12)^2)*as.vector(xx_R1%*%Beta_1[,6])*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(y11*Conj(y12)*(a[1]) + Conj(y11*Conj(y12)*(a[1])))*t(xx_R1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(-ye2*Conj(ye3) - ye3*Conj(ye2))*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(ye2)^2)*as.vector(xx_Re%*%Beta_1[,6])*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(ye1*Conj(ye2)*(a[eee]) + Conj(ye1*Conj(ye2)*(a[eee])))*t(xx_Re))
                        
                        temp_mat_71 <- xx_I*ik7
                        temp_mat_72 <- ck7*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,1])))
                        gr7 <- gr7 + colSums((temp_mat_71 + temp_mat_72)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(Im(-y11*Conj(y12) + y12*Conj(y11)))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2)*as.vector(xx_I1%*%Beta_2[,1])*t(xx_I1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(Im(-ye1*Conj(ye2) + ye2*Conj(ye1)))*t(xx_Ie) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(ye1)^2)*as.vector(xx_Ie%*%Beta_2[,1])*t(xx_Ie))
                        temp_mat_81 <- xx_I*ik8
                        temp_mat_82 <- ck8*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,2])))
                        temp_mat_83 <- xx_I*dk8
                        gr8 <- gr8 + colSums((temp_mat_81 + temp_mat_82 + temp_mat_83)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Im(-y11*Conj(y13) + y13*Conj(y11)))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2)*as.vector(xx_I1%*%Beta_2[,2])*t(xx_I1) +
                                              (1i)*as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y11*Conj(y12)*(b[1]) + Conj(y11*Conj(y12)*(b[1])))*t(xx_I1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(Im(-ye1*Conj(ye3) + ye3*Conj(ye1)))*t(xx_Ie) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(ye1)^2)*as.vector(xx_Ie%*%Beta_2[,2])*t(xx_Ie) +
                                              (1i)*as.vector(exp(-xx_Re%*%Beta_1[,3]))*(-ye1*Conj(ye2)*(b[eee]) + Conj(ye1*Conj(ye2)*(b[eee])))*t(xx_Ie))
                        
                        temp_mat_91 <- xx_I*ik9
                        temp_mat_92 <- ck9*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,3])))
                        temp_mat_93 <- xx_I*dk9
                        gr9 <- gr9 + colSums((temp_mat_91 + temp_mat_92 + temp_mat_93)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Im(-y12*Conj(y13) + y13*Conj(y12)))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y12)^2)*as.vector(xx_I1%*%Beta_2[,3])*t(xx_I1) +
                                              (1i)*as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y11*Conj(y12)*(a[1]) + Conj(y11*Conj(y12)*(a[1])))*t(xx_I1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(Im(-ye2*Conj(ye3) + ye3*Conj(ye2)))*t(xx_Ie) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(ye2)^2)*as.vector(xx_Ie%*%Beta_2[,3])*t(xx_Ie) +
                                              (1i)*as.vector(exp(-xx_Re%*%Beta_1[,3]))*(-ye1*Conj(ye2)*(a[eee]) + Conj(ye1*Conj(ye2)*(a[eee])))*t(xx_Ie))
                        #==============
                        # Hession
                        #==============
                        bigmat_h11 <- kronecker(as.vector(abs(y1)^2*exp(-crossprod(t(xx_R),Beta_1[,1])))*xx_R, t(rep(1,nBeta)))
                        coefmat_h11 <- repmat(xx_R,1,nBeta)
                        h11 <- h11 + matrix(colSums(bigmat_h11*coefmat_h11),nBeta,nBeta) +
                                0.5*(abs(y11)^2*as.vector(exp(-xx_R1%*%Beta_1[,1]))*xx_R1%*%t(xx_R1)) +
                                0.5*(abs(ye1)^2*as.vector(exp(-xx_Re%*%Beta_1[,1]))*xx_Re%*%t(xx_Re))
                        
                        bigmat_h22 <- kronecker(as.vector(abs(y2-theta[1,2:nfreq]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))
                                                *xx_R,t(rep(1,nBeta)))
                        coefmat_h22 <- repmat(xx_R,1,nBeta)
                        h22 <- h22 + matrix(colSums(bigmat_h22*coefmat_h22),nBeta,nBeta) +
                                0.5*(abs(y12-theta[1,1]*y11)^2*as.vector(exp(-xx_R1%*%Beta_1[,2]))*xx_R1%*%t(xx_R1)) +
                                0.5*(abs(ye2-theta[1,eee]*ye1)^2*as.vector(exp(-xx_Re%*%Beta_1[,2]))*xx_Re%*%t(xx_Re))
                        
                        bigmat_h33 <- kronecker(as.vector(abs(y3-theta[2,2:nfreq]*y1 - theta[3,2:nfreq]*y2)^2*exp(-crossprod(t(xx_R),Beta_1[,3])))
                                                *xx_R,t(rep(1,nBeta)))
                        coefmat_h33 <- repmat(xx_R,1,nBeta)
                        h33 <- h33 + matrix(colSums(bigmat_h33*coefmat_h33),nBeta,nBeta) +       
                                0.5*(abs(y13-theta[2,1]*y11 - theta[3,1]*y12)^2*as.vector(exp(-xx_R1%*%Beta_1[,3]))*xx_R1%*%t(xx_R1)) +
                                0.5*(abs(ye3-theta[2,eee]*ye1 - theta[3,eee]*ye2)^2*as.vector(exp(-xx_Re%*%Beta_1[,3]))*xx_Re%*%t(xx_Re))
                        
                        bigmat_h44 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,2]))*ck4)*xx_R, t(rep(1,nBeta)))
                        coefmat_h44 <- repmat(xx_R,1,nBeta)
                        h44 <- h44 + matrix(colSums(bigmat_h44*coefmat_h44),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2])*2*abs(y11)^2)*xx_R1%*%t(xx_R1)) +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,2])*2*abs(ye1)^2)*xx_Re%*%t(xx_Re))
                        
                        bigmat_h55 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,3]))*ck5)*xx_R, t(rep(1,nBeta)))
                        coefmat_h55 <- repmat(xx_R,1,nBeta)
                        h55 <- h55 + matrix(colSums(bigmat_h55*coefmat_h55),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3])*2*abs(y11)^2)*xx_R1%*%t(xx_R1)) +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3])*2*abs(ye1)^2)*xx_Re%*%t(xx_Re))
                        
                        bigmat_h66 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,3]))*ck6)*xx_R, t(rep(1,nBeta)))
                        coefmat_h66 <- repmat(xx_R,1,nBeta)
                        h66 <- h66 + matrix(colSums(bigmat_h66*coefmat_h66),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3])*2*abs(y12)^2)*xx_R1%*%t(xx_R1))  +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3])*2*abs(ye2)^2)*xx_Re%*%t(xx_Re)) 
                        
                        bigmat_h77 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,2]))*ck7)*xx_I, t(rep(1,nBeta)))
                        coefmat_h77 <- repmat(xx_I,1,nBeta)
                        h77 <- h77 + matrix(colSums(bigmat_h77*coefmat_h77),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2])*2*abs(y11)^2)*xx_I1%*%t(xx_I1)) +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,2])*2*abs(ye1)^2)*xx_Ie%*%t(xx_Ie))
                        
                        bigmat_h88 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,3]))*ck8)*xx_I, t(rep(1,nBeta)))
                        coefmat_h88 <- repmat(xx_I,1,nBeta)
                        h88 <- h88 + matrix(colSums(bigmat_h88*coefmat_h88),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3])*2*abs(y11)^2)*xx_I1%*%t(xx_I1))  +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3])*2*abs(ye1)^2)*xx_Ie%*%t(xx_Ie))
                        
                        bigmat_h99 <- kronecker(as.vector(exp(-crossprod(t(xx_R), Beta_1[,3]))*ck9)*xx_I, t(rep(1,nBeta)))
                        coefmat_h99 <- repmat(xx_I,1,nBeta)
                        h99 <- h99 + matrix(colSums(bigmat_h99*coefmat_h99),nBeta,nBeta) +
                                0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3])*2*abs(y12)^2)*xx_I1%*%t(xx_I1)) +
                                0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3])*2*abs(ye2)^2)*xx_Ie%*%t(xx_Ie)) 
                        
                        bigmat_h42_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck4*(crossprod(t(xx_R),Beta_1[,4])))
                                                  *xx_R, t(rep(1,nBeta)))
                        coefmat_h42_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h42_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*rk4)*xx_R, t(rep(1,nBeta)))
                        coefmat_h42_2 <- repmat(xx_R, 1,nBeta)
                        h42 <- h42 + t(matrix(colSums(bigmat_h42_1*coefmat_h42_1 + 
                                        bigmat_h42_2*coefmat_h42_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2*as.vector(xx_R1%*%Beta_1[,4])*xx_R1 +
                                        (-y11*Conj(y12) - y12*Conj(y11))*xx_R1)%*%t(xx_R1)) +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(ye1)^2*as.vector(xx_Re%*%Beta_1[,4])*xx_Re +
                                        (-ye1*Conj(ye2) - ye2*Conj(ye1))*xx_Re)%*%t(xx_Re))
                        
                        bigmat_h53_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ck5*(crossprod(t(xx_R),Beta_1[,5])))
                                                  *xx_R, t(rep(1,nBeta)))
                        coefmat_h53_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h53_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*rk5)*xx_R, t(rep(1,nBeta)))
                        coefmat_h53_2 <- repmat(xx_R, 1,nBeta)
                        bigmat_h53_3 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*dk5)*xx_R, t(rep(1,nBeta)))
                        coefmat_h53_3 <- repmat(xx_R, 1,nBeta)                        
                        h53 <- h53 + t(matrix(colSums(bigmat_h53_1*coefmat_h53_1 + 
                                        bigmat_h53_2*coefmat_h53_2 +
                                        bigmat_h53_3*coefmat_h53_3),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2*as.vector(xx_R1%*%Beta_1[,5])*xx_R1 +
                                        (-y11*Conj(y13) - y13*Conj(y11))*xx_R1 + 
                                        (y12*Conj(y11)*b[1] + Conj(y12*Conj(y11)*b[1]))*xx_R1)%*%t(xx_R1)) +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(ye1)^2*as.vector(xx_Re%*%Beta_1[,5])*xx_Re +
                                        (-ye1*Conj(ye3) - ye3*Conj(ye1))*xx_Re + 
                                        (ye2*Conj(ye1)*b[eee] + Conj(ye2*Conj(ye1)*b[eee]))*xx_Re)%*%t(xx_Re))
                        
                        bigmat_h56_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y1)*y2)*xx_R, t(rep(1,nBeta)))
                        coefmat_h56_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h56_2 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y2)*y1)*xx_R, t(rep(1,nBeta)))
                        coefmat_h56_2 <- repmat(xx_R, 1,nBeta)
                        h56 <- Re(h56 + t(matrix(colSums(bigmat_h56_1*coefmat_h56_1 + 
                                        bigmat_h56_2*coefmat_h56_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Conj(y11)*y12*xx_R1 +
                                        Conj(y12)*y11*xx_R1)%*%t(xx_R1)) +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(Conj(ye1)*ye2*xx_Re +
                                        Conj(ye2)*ye1*xx_Re)%*%t(xx_Re))) 
                        
                        bigmat_h59_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y1)*y2)*xx_R, t(rep(1,nBeta)))
                        coefmat_h59_1 <- repmat(xx_I,1,nBeta)
                        bigmat_h59_2 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y2)*y1)*xx_R, t(rep(1,nBeta)))
                        coefmat_h59_2 <- repmat(xx_I, 1,nBeta)
                        h59 <- Im(h59 - t(matrix(colSums(bigmat_h59_1*coefmat_h59_1 - 
                                        bigmat_h59_2*coefmat_h59_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Conj(y11)*y12*xx_R1 +
                                        Conj(y12)*y11*xx_R1)%*%t(xx_I1)) +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(Conj(ye1)*ye2*xx_Re +
                                        Conj(ye2)*ye1*xx_Re)%*%t(xx_Ie)))  
                        
                        bigmat_h63_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ck6*(crossprod(t(xx_R),Beta_1[,6])))
                                                  *xx_R, t(rep(1,nBeta)))
                        coefmat_h63_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h63_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*rk6)*xx_R, t(rep(1,nBeta)))
                        coefmat_h63_2 <- repmat(xx_R, 1,nBeta)
                        bigmat_h63_3 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*dk6)*xx_R, t(rep(1,nBeta)))
                        coefmat_h63_3 <- repmat(xx_R, 1,nBeta)                        
                        h63 <- h63 + t(matrix(colSums(bigmat_h63_1*coefmat_h63_1 + 
                                        bigmat_h63_2*coefmat_h63_2 +
                                        bigmat_h63_3*coefmat_h63_3),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2*as.vector(xx_R1%*%Beta_1[,6])*xx_R1 +
                                        (-y12*Conj(y13) - y13*Conj(y12))*xx_R1 + 
                                        (y12*Conj(y11)*a[1] + Conj(y12*Conj(y11)*a[1]))*xx_R1)%*%t(xx_R1))  +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(ye1)^2*as.vector(xx_Re%*%Beta_1[,6])*xx_Re +
                                        (-ye2*Conj(ye3) - ye3*Conj(ye2))*xx_Re + 
                                        (ye2*Conj(ye1)*a[eee] + Conj(ye2*Conj(ye1)*a[eee]))*xx_Re)%*%t(xx_Re)) 
                                        
                        bigmat_h68_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y2)*y1)*xx_R, t(rep(1,nBeta)))
                        coefmat_h68_1 <- repmat(xx_I,1,nBeta)
                        bigmat_h68_2 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y1)*y2)*xx_R, t(rep(1,nBeta)))
                        coefmat_h68_2 <- repmat(xx_I, 1,nBeta)
                        h68 <- Im(h68 - t(matrix(colSums(bigmat_h68_1*coefmat_h68_1 -
                                        bigmat_h68_2*coefmat_h68_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Conj(y12)*y11*xx_R1 +
                                        Conj(y11)*y12*xx_R1)%*%t(xx_I1))  + 
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(Conj(ye2)*ye1*xx_Re +
                                        Conj(ye1)*ye2*xx_Re)%*%t(xx_Ie)))
                        
                        bigmat_h72_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ck7*(crossprod(t(xx_I),Beta_2[,1])))
                                                  *xx_I, t(rep(1,nBeta)))
                        coefmat_h72_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h72_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,2]))*ik7)*xx_I, t(rep(1,nBeta)))  
                        coefmat_h72_2 <- repmat(xx_R, 1,nBeta)
                        h72 <- h72 + t(matrix(colSums(bigmat_h72_1*coefmat_h72_1 + 
                                        bigmat_h72_2*coefmat_h72_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y11)^2*as.vector(xx_I1%*%Beta_2[,1])*xx_I1 +
                                        Im(-y11*Conj(y12) + y12*Conj(y11))*xx_I1)%*%t(xx_R1)) +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(ye1)^2*as.vector(xx_Ie%*%Beta_2[,1])*xx_Ie +
                                        Im(-ye1*Conj(ye2) + ye2*Conj(ye1))*xx_Ie)%*%t(xx_Re))
                        
                        bigmat_h83_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ck8*(crossprod(t(xx_I),Beta_2[,2])))
                                                  *xx_I, t(rep(1,nBeta)))
                        coefmat_h83_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h83_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ik8)*xx_I, t(rep(1,nBeta)))
                        coefmat_h83_2 <- repmat(xx_R, 1,nBeta)
                        bigmat_h83_3 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*dk8)*xx_I, t(rep(1,nBeta)))
                        coefmat_h83_3 <- repmat(xx_R, 1,nBeta)                        
                        h83 <- h83 + t(matrix(colSums(bigmat_h83_1*coefmat_h83_1 + 
                                        bigmat_h83_2*coefmat_h83_2 +
                                        bigmat_h83_3*coefmat_h83_3),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y11)^2*as.vector(xx_I1%*%Beta_2[,2])*xx_I1 +
                                        (-y11*Conj(y13) + y13*Conj(y11))*xx_I1 + 
                                        Im(-y12*Conj(y11)*b[1] + Conj(y12*Conj(y11)*b[1]))*xx_I1)%*%t(xx_R1))  +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(ye1)^2*as.vector(xx_Ie%*%Beta_2[,2])*xx_Ie +
                                        (-ye1*Conj(ye3) + ye3*Conj(ye1))*xx_Ie + 
                                        Im(-ye2*Conj(ye1)*b[eee] + Conj(ye2*Conj(ye1)*b[eee]))*xx_Ie)%*%t(xx_Re))           
                                        
                        bigmat_h93_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ck9*(crossprod(t(xx_I),Beta_2[,3])))
                                                  *xx_I, t(rep(1,nBeta)))
                        coefmat_h93_1 <- repmat(xx_R,1,nBeta)
                        bigmat_h93_2 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*ik9)*xx_I, t(rep(1,nBeta)))
                        coefmat_h93_2 <- repmat(xx_R, 1,nBeta)
                        bigmat_h93_3 <- kronecker( as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*dk9)*xx_I, t(rep(1,nBeta)))
                        coefmat_h93_3 <- repmat(xx_R, 1,nBeta)                        
                        h93 <- h93 + t(matrix(colSums(bigmat_h93_1*coefmat_h93_1 + 
                                        bigmat_h93_2*coefmat_h93_2 +
                                        bigmat_h93_3*coefmat_h93_3),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y12)^2*as.vector(xx_I1%*%Beta_2[,3])*xx_I1 +
                                        (-y12*Conj(y13) + y13*Conj(y12))*xx_I1 + 
                                        Im(-y11*Conj(y12)*a[1] + Conj(y11*Conj(y12)*a[1]))*xx_I1)%*%t(xx_R1)) +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(ye2)^2*as.vector(xx_Ie%*%Beta_2[,3])*xx_Ie +
                                        (-ye2*Conj(ye3) + ye3*Conj(ye2))*xx_Ie + 
                                        Im(-ye1*Conj(ye2)*a[eee] + Conj(ye1*Conj(ye2)*a[eee]))*xx_Ie)%*%t(xx_Re))   
                                        
                        bigmat_h98_1 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y1)*y2)*xx_I, t(rep(1,nBeta)))
                        coefmat_h98_1 <- repmat(xx_I,1,nBeta)
                        bigmat_h98_2 <- kronecker(as.vector(exp(-crossprod(t(xx_R),Beta_1[,3]))*Conj(y2)*y1)*xx_I, t(rep(1,nBeta)))
                        coefmat_h98_2 <- repmat(xx_I, 1,nBeta)
                        h98 <- Re(h98 + t(matrix(colSums(bigmat_h98_1*coefmat_h98_1 + 
                                        bigmat_h98_2*coefmat_h98_2),nBeta,nBeta)) +
                                        0.5*(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Conj(y11)*y12*xx_I1 +
                                        Conj(y12)*y11*xx_I1)%*%t(xx_I1))  +
                                        0.5*(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(Conj(ye1)*ye2*xx_Ie +
                                        Conj(ye2)*ye1*xx_Ie)%*%t(xx_Ie)))  
                        
                        h24=t(h42); h35=t(h53); h65=t(h56); h95=t(h59); h36=t(h63);
                        h86=t(h68); h27=t(h72); h38=t(h83); h39=t(h93); h89=t(h98);
                }
                ze <- matrix(0,nBeta,nBeta)
                h1 <- cbind(h11, repmat(ze,1,8))
                h2 <- cbind(ze, h22, ze, -h24, ze, ze, -h27, ze, ze)
                h3 <- cbind(ze, ze, h33, ze, -h35, -h36, ze, -h38, -h39)
                h4 <- cbind(ze, -h42, ze, h44, repmat(ze,1,5))
                h5 = cbind(ze,ze,-h53,ze,h55,h56,ze,ze,h59)
                h6 = cbind(ze,ze,-h63,ze,h65,h66,ze,h68,ze)
                h7 = cbind(ze,-h72,ze,ze,ze,ze,h77,ze,ze)
                h8 = cbind(ze,ze,-h83,ze,ze,h86,ze,h88,h89)
                h9 = cbind(ze,ze,-h93,ze,h95,ze,ze,h98,h99)
                
                gr <- c(gr1,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9); h <- rbind(h1,h2,h3,h4,h5,h6,h7,h8,h9); f <- -Re(f)
                gr_index <- (1:(dimen^2*nBeta))*kronecker(chol_index[Phi_temp,],rep(1,nBeta))
                gr_index = gr_index[which(gr_index!=0)]
                gr = Re(gr[gr_index]); h = Re(h[gr_index,gr_index])
        }
        list(value = f, gradient = gr, hessian = h)
}
