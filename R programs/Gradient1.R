
Gradient1 <- function(yobs_tmp, chol_index, Phi_temp, tau_temp,
                         Beta_temp, sigmasqalpha, nbasis, global_info){
        
        dimen <- global_info$dimen
        
        nBeta <- nbasis + 1
        Beta_1 <- matrix(0,nBeta,(dimen + dimen*(dimen-1)/2))
        Beta_2 <- matrix(0,nBeta,dimen*(dimen-1)/2)
        Beta_1[,] <- Beta_temp[,1:(dimen + dimen*(dimen-1)/2)]
        Beta_2[,] <- Beta_temp[,-(1:(dimen + dimen*(dimen-1)/2))]
        
        n <- dim(yobs_tmp)[1]
        nfreq <- floor(n/2); tt <- (0:nfreq)/(2*nfreq)
        yy <- mvfft(yobs_tmp)/sqrt(n); y <- yy[1:(nfreq+1),]; nf <- dim(y)[1]; eee <- dim(y)[1]
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
        
        if (dimen==2){
                
                gr1 <- matrix(0,nBeta,1);  gr2 <- matrix(0,nBeta,1); gr3 <- matrix(0,nBeta,1); gr4 <- matrix(0,nBeta,1)
                gr1[1] <- Beta_1[1,1]/sigmasqalpha; gr1[2:nBeta,1] <- Beta_1[2:nBeta,1]/tau_temp[1]
                gr2[1] <- Beta_1[1,2]/sigmasqalpha; gr1[2:nBeta,1] <- Beta_1[2:nBeta,2]/tau_temp[2]
                gr3[1] <- Beta_1[1,3]/sigmasqalpha; gr1[2:nBeta,1] <- Beta_1[2:nBeta,3]/tau_temp[3]
                gr4[1:nBeta,1] <- Beta_2[1:nBeta,1]/tau_temp[4]
                
                if(n%%2==1){
                        #==============
                        # gradient
                        #==============
                        rk <- -y[2:eee,1]*Conj(y[2:eee,2]) - y[2:eee,2]*Conj(y[2:eee,1])
                        ik <- 1i*(-y[2:eee,1]*Conj(y[2:eee,2]) + y[2:eee,2]*Conj(y[2:eee,1]))
                        ck <- 2*abs(y[2:eee,1])^2
                        
                        gr1 <- gr1 + crossprod(xx_r[2:eee,], (1-abs(y[2:eee,1])^2*exp(-xx_r[2:eee,]%*%Beta_1[,1]))) +
                                0.5*(xx_r[1,]*(1-abs(y[1,1])^2*exp(-xx_r[1,]%*%Beta_1[,1])))
                        
                        gr2 <- gr2 + crossprod(xx_r[2:eee,], (1-abs(y[2:eee,2] - theta[2:eee]*y[2:eee,1])^2*exp(-xx_r[2:eee,]%*%Beta_1[,2]))) +
                                0.5*(xx_r[1,]*(1-abs(y[1,2] - t(theta[1])*y[1,1])^2*exp(-xx_r[1,]%*%Beta_1[,2])))
                        
                        temp_mat_31 <- xx_r[2:eee,]*rk
                        temp_mat_32 <- ck*(xx_r[2:eee,]*as.vector(crossprod(t(xx_r[2:eee,]),Beta_1[,3])))
                        gr3 <- gr3 + colSums((temp_mat_31 + temp_mat_32)*as.vector(exp(-crossprod(t(xx_r[2:eee,]),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_r[1,]%*%Beta_1[,2]))*(-y[1,1]*Conj(y[1,2]) - y[1,2]*Conj(y[1,1]))*t(xx_r[1,]) +
                                              as.vector(exp(-xx_r[1,]%*%Beta_1[,2]))*(2*abs(y[1,1])^2)*as.vector(xx_r[1,]%*%Beta_1[,3])*t(xx_r[1,]))
                        
                        temp_mat_41 <- ik*xx_i[2:eee,]
                        temp_mat_42 <- ck*(xx_i[2:eee,]*as.vector(crossprod(t(xx_i[2:eee,]),Beta_2[,1])))
                        gr4 <- gr4 + colSums((temp_mat_41 + temp_mat_42)*as.vector(exp(-crossprod(t(xx_r[2:eee,]),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_r[1,]%*%Beta_1[,2]))*(1i*(-y[1,1]*Conj(y[1,2]) + y[1,2]*Conj(y[1,1])))*t(xx_i[1,]) +
                                              as.vector(exp(-xx_r[1,]%*%Beta_1[,2]))*(2*abs(y[1,1])^2)*as.vector(xx_i[1,]%*%Beta_2[,1])*t(xx_i[1,]))            

                }else{
                        #==============
                        # gradient
                        #==============
                        rk <- -y[2:nfreq,1]*Conj(y[2:nfreq,2]) - y[2:nfreq,2]*Conj(y[2:nfreq,1])
                        ik <- 1i*(-y[2:nfreq,1]*Conj(y[2:nfreq,2]) + y[2:nfreq,2]*Conj(y[2:nfreq,1]))
                        ck <- 2*abs(y[2:nfreq,1])^2
                        
                        gr1 <- gr1 + crossprod(xx_r[2:nfreq,], (1-abs(y[2:nfreq,1])^2*exp(-xx_r[2:nfreq,]%*%Beta_1[,1]))) +
                                0.5*(xx_r[1,]*(1-abs(y[1,1])^2*exp(-xx_r[1,]%*%Beta_1[,1]))) +
                                0.5*(xx_r[eee,]*(1-abs(y[eee,1])^2*exp(-xx_r[eee,]%*%Beta_1[,1])))
                        
                        gr2 <- gr2 + crossprod(xx_r[2:nfreq,], (1-abs(y[2:nfreq,2] - theta[2:nfreq]*y[2:nfreq,1])^2*exp(-xx_r[2:nfreq,]%*%Beta_1[,2]))) +
                                0.5*(xx_r[1,]*(1-abs(y[1,2] - t(theta[1])*y[1,1])^2*exp(-xx_r[1,]%*%Beta_1[,2]))) +
                                0.5*(xx_r[eee,]*(1-abs(y[eee,2] - t(theta[eee])*y[eee,1])^2*exp(-xx_r[eee,]%*%Beta_1[,2])))
                        
                        temp_mat_31 <- xx_r[2:nfreq,]*rk
                        temp_mat_32 <- ck*(xx_r[2:nfreq,]*as.vector(crossprod(t(xx_r[2:nfreq,]),Beta_1[,3])))
                        gr3 <- gr3 + colSums((temp_mat_31 + temp_mat_32)*as.vector(exp(-crossprod(t(xx_r[2:nfreq,]),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_r[1,]%*%Beta_1[,2]))*(-y[1,1]*Conj(y[1,2]) - y[1,2]*Conj(y[1,1]))*t(xx_r[1,]) +
                                              as.vector(exp(-xx_r[1,]%*%Beta_1[,2]))*(2*abs(y[1,1])^2)*as.vector(xx_r[1,]%*%Beta_1[,3])*t(xx_r[1,])) +
                                0.5*t(as.vector(exp(-xx_r[eee,]%*%Beta_1[,2]))*(-y[eee,1]*Conj(y[eee,2]) - y[eee,2]*Conj(y[eee,1]))*t(xx_r[eee,]) +
                                              as.vector(exp(-xx_r[eee,]%*%Beta_1[,2]))*(2*abs(y[eee,1])^2)*as.vector(xx_r[eee,]%*%Beta_1[,3])*t(xx_r[eee,]))
                        
                        temp_mat_41 <- ik*xx_i[2:nfreq,]
                        temp_mat_42 <- ck*(xx_i[2:nfreq,]*as.vector(crossprod(t(xx_i[2:nfreq,]),Beta_2[,1])))
                        gr4 <- gr4 + colSums((temp_mat_41 + temp_mat_42)*as.vector(exp(-crossprod(t(xx_r[2:nfreq,]),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_r[1,]%*%Beta_1[,2]))*(1i*(-y[1,1]*Conj(y[1,2]) + y[1,2]*Conj(y[1,1])))*t(xx_i[1,]) +
                                              as.vector(exp(-xx_r[1,]%*%Beta_1[,2]))*(2*abs(y[1,1])^2)*as.vector(xx_i[1,]%*%Beta_2[,1])*t(xx_i[1,])) +
                                0.5*t(as.vector(exp(-xx_r[eee,]%*%Beta_1[,2]))*(1i*(-y[eee,1]*Conj(y[eee,2]) + y[eee,2]*Conj(y[eee,1])))*t(xx_i[eee,]) +
                                              as.vector(exp(-xx_r[eee,]%*%Beta_1[,2]))*(2*abs(y[eee,1])^2)*as.vector(xx_i[eee,]%*%Beta_2[,1])*t(xx_i[eee,]))
                }
                gr <- c(gr1,gr2,gr3,gr4) 
                gr_index <- (1:(4*nBeta))* kron(chol_index[Phi_temp,1:4],rep(1,nBeta))
                gr_index <- gr_index[which(gr_index!=0)]
                gr <- Re(gr[gr_index])
                gr <- -gr
        } else if (dimen==3){
                
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
                
                if(n%%2==1){
                        y1 <- y[2:eee,1];y2 <- y[2:eee,2]; y3 <- y[2:eee,3]
                        xx_R <- xx_r[2:eee,];xx_I <- xx_i[2:eee,]
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
                                0.5*(xx_R1*(1-abs(y[1,1])^2*exp(-xx_R1%*%Beta_1[,1])))
                        gr2 <- gr2 + crossprod(xx_R, (1-abs(y2 - theta[1,2:eee]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*(xx_R1*(1-abs(y[1,2] - t(theta[1,1])*y[1,1])^2*exp(-xx_R1%*%Beta_1[,2])))
                        gr3 <- gr3 + crossprod(xx_R, (1-abs(y3 - theta[2,2:eee]*y1 - theta[3,2:eee]*y2 )^2*exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*(xx_R1*(1-abs(y[1,3] - t(theta[2,1])*y[1,1] - t(theta[3,1])*y[1,2])^2*exp(-xx_R1%*%Beta_1[,3])))
                        temp_mat_41 <- xx_R*rk4
                        temp_mat_42 <- ck4*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,4])))
                        gr4 <- gr4 + colSums((temp_mat_41 + temp_mat_42)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(-y[1,1]*Conj(y[1,2]) - y[1,2]*Conj(y[1,1]))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y[1,1])^2)*as.vector(xx_R1%*%Beta_1[,4])*t(xx_R1))
                        temp_mat_51 <- xx_R*rk5
                        temp_mat_52 <- ck5*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,5])))
                        temp_mat_53 <- xx_R*dk5
                        gr5 <- gr5 + colSums((temp_mat_51 + temp_mat_52 + temp_mat_53)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y[1,1]*Conj(y[1,3]) - y[1,3]*Conj(y[1,1]))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y[1,1])^2)*as.vector(xx_R1%*%Beta_1[,5])*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(y[1,2]*Conj(y[1,1])*(b[1]) + Conj(y[1,2]*Conj(y[1,1])*(b[1])))*t(xx_R1))
                        temp_mat_61 <- xx_R*rk6
                        temp_mat_62 <- ck6*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,6])))
                        temp_mat_63 <- xx_R*dk6
                        gr6 <- gr6 + colSums((temp_mat_61 + temp_mat_62 + temp_mat_63)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y[1,2]*Conj(y[1,3]) - y[1,3]*Conj(y[1,2]))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y[1,2])^2)*as.vector(xx_R1%*%Beta_1[,6])*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(y[1,1]*Conj(y[1,2])*(a[1]) + Conj(y[1,1]*Conj(y[1,2])*(a[1])))*t(xx_R1))
                        temp_mat_71 <- xx_I*ik7
                        temp_mat_72 <- ck7*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,1])))
                        gr7 <- gr7 + colSums((temp_mat_71 + temp_mat_72)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(Im(-y[1,1]*Conj(y[1,2]) + y[1,2]*Conj(y[1,1])))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y[1,1])^2)*as.vector(xx_I1%*%Beta_2[,1])*t(xx_I1))
                        temp_mat_81 <- xx_I*ik8
                        temp_mat_82 <- ck8*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,2])))
                        temp_mat_83 <- xx_I*dk8
                        gr8 <- gr8 + colSums((temp_mat_81 + temp_mat_82 + temp_mat_83)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Im(-y[1,1]*Conj(y[1,3]) + y[1,3]*Conj(y[1,1])))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y[1,1])^2)*as.vector(xx_I1%*%Beta_2[,2])*t(xx_I1) +
                                              (1i)*as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y[1,1]*Conj(y[1,2])*(b[1]) + Conj(y[1,1]*Conj(y[1,2])*(b[1])))*t(xx_I1))                        
                        temp_mat_91 <- xx_I*ik9
                        temp_mat_92 <- ck9*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,3])))
                        temp_mat_93 <- xx_I*dk9
                        gr9 <- gr9 + colSums((temp_mat_91 + temp_mat_92 + temp_mat_93)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Im(-y[1,2]*Conj(y[1,3]) + y[1,3]*Conj(y[1,2])))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y[1,2])^2)*as.vector(xx_I1%*%Beta_2[,3])*t(xx_I1) +
                                              (1i)*as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y[1,1]*Conj(y[1,2])*(a[1]) + Conj(y[1,1]*Conj(y[1,2])*(a[1])))*t(xx_I1)) 
      
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
                                0.5*(xx_R1*(1-abs(y[1,1])^2*exp(-xx_R1%*%Beta_1[,1])))+
                                0.5*(xx_Re*(1-abs(y[eee,1])^2*exp(-xx_Re%*%Beta_1[,1])))
                        gr2 <- gr2 + crossprod(xx_R, (1-abs(y2 - theta[1,2:nfreq]*y1)^2*exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*(xx_R1*(1-abs(y[1,2] - theta[1,1]*y[1,1])^2*exp(-xx_R1%*%Beta_1[,2])))+
                                0.5*(xx_Re*(1-abs(y[eee,2] - theta[1,eee]*y[eee,1])^2*exp(-xx_Re%*%Beta_1[,2])))
                        gr3 <- gr3 + crossprod(xx_R, (1-abs(y3 - theta[2,2:nfreq]*y1 - theta[3,2:nfreq]*y2 )^2*exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*(xx_R1*(1-abs(y[1,3] - t(theta[2,1])*y[1,1] - t(theta[3,1])*y[1,2])^2*exp(-xx_R1%*%Beta_1[,3]))) +
                                0.5*(xx_Re*(1-abs(y[eee,3] - t(theta[2,eee])*y[eee,1] - t(theta[3,eee])*y[eee,2])^2*exp(-xx_Re%*%Beta_1[,3])))
                        temp_mat_41 <- xx_R*rk4
                        temp_mat_42 <- ck4*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,4])))
                        gr4 <- gr4 + colSums((temp_mat_41 + temp_mat_42)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(-y[1,1]*Conj(y[1,2]) - y[1,2]*Conj(y[1,1]))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y[1,1])^2)*as.vector(xx_R1%*%Beta_1[,4])*t(xx_R1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(-y[eee,1]*Conj(y[eee,2]) - y[eee,2]*Conj(y[eee,1]))*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(y[eee,1])^2)*as.vector(xx_Re%*%Beta_1[,4])*t(xx_Re))
                        temp_mat_51 <- xx_R*rk5
                        temp_mat_52 <- ck5*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,5])))
                        temp_mat_53 <- xx_R*dk5
                        gr5 <- gr5 + colSums((temp_mat_51 + temp_mat_52 + temp_mat_53)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y[1,1]*Conj(y[1,3]) - y[1,3]*Conj(y[1,1]))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y[1,1])^2)*as.vector(xx_R1%*%Beta_1[,5])*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(y[1,2]*Conj(y[1,1])*(b[1]) + Conj(y[1,2]*Conj(y[1,1])*(b[1])))*t(xx_R1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(-y[eee,1]*Conj(y[eee,3]) - y[eee,3]*Conj(y[eee,1]))*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(y[eee,1])^2)*as.vector(xx_Re%*%Beta_1[,5])*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(y[eee,2]*Conj(y[eee,1])*(b[eee]) + Conj(y[eee,2]*Conj(y[eee,1])*(b[eee])))*t(xx_Re))
                        temp_mat_61 <- xx_R*rk6
                        temp_mat_62 <- ck6*(xx_R*as.vector(crossprod(t(xx_R),Beta_1[,6])))
                        temp_mat_63 <- xx_R*dk6
                        gr6 <- gr6 + colSums((temp_mat_61 + temp_mat_62 + temp_mat_63)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y[1,2]*Conj(y[1,3]) - y[1,3]*Conj(y[1,2]))*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y[1,2])^2)*as.vector(xx_R1%*%Beta_1[,6])*t(xx_R1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(y[1,1]*Conj(y[1,2])*(a[1]) + Conj(y[1,1]*Conj(y[1,2])*(a[1])))*t(xx_R1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(-y[eee,2]*Conj(y[eee,3]) - y[eee,3]*Conj(y[eee,2]))*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(y[eee,2])^2)*as.vector(xx_Re%*%Beta_1[,6])*t(xx_Re) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(y[eee,1]*Conj(y[eee,2])*(a[eee]) + Conj(y[eee,1]*Conj(y[eee,2])*(a[eee])))*t(xx_Re))
                        
                        temp_mat_71 <- xx_I*ik7
                        temp_mat_72 <- ck7*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,1])))
                        gr7 <- gr7 + colSums((temp_mat_71 + temp_mat_72)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,2])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,2]))*(Im(-y[1,1]*Conj(y[1,2]) + y[1,2]*Conj(y[1,1])))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,2]))*(2*abs(y[1,1])^2)*as.vector(xx_I1%*%Beta_2[,1])*t(xx_I1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,2]))*(Im(-y[eee,1]*Conj(y[eee,2]) + y[eee,2]*Conj(y[eee,1])))*t(xx_Ie) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,2]))*(2*abs(y[eee,1])^2)*as.vector(xx_Ie%*%Beta_2[,1])*t(xx_Ie))
                        temp_mat_81 <- xx_I*ik8
                        temp_mat_82 <- ck8*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,2])))
                        temp_mat_83 <- xx_I*dk8
                        gr8 <- gr8 + colSums((temp_mat_81 + temp_mat_82 + temp_mat_83)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Im(-y[1,1]*Conj(y[1,3]) + y[1,3]*Conj(y[1,1])))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y[1,1])^2)*as.vector(xx_I1%*%Beta_2[,2])*t(xx_I1) +
                                              (1i)*as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y[1,1]*Conj(y[1,2])*(b[1]) + Conj(y[1,1]*Conj(y[1,2])*(b[1])))*t(xx_I1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(Im(-y[eee,1]*Conj(y[eee,3]) + y[eee,3]*Conj(y[eee,1])))*t(xx_Ie) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(y[eee,1])^2)*as.vector(xx_Ie%*%Beta_2[,2])*t(xx_Ie) +
                                              (1i)*as.vector(exp(-xx_Re%*%Beta_1[,3]))*(-y[eee,1]*Conj(y[eee,2])*(b[eee]) + Conj(y[eee,1]*Conj(y[eee,2])*(b[eee])))*t(xx_Ie))
                        
                        temp_mat_91 <- xx_I*ik9
                        temp_mat_92 <- ck9*(xx_I*as.vector(crossprod(t(xx_I),Beta_2[,3])))
                        temp_mat_93 <- xx_I*dk9
                        gr9 <- gr9 + colSums((temp_mat_91 + temp_mat_92 + temp_mat_93)*as.vector(exp(-crossprod(t(xx_R),Beta_1[,3])))) +
                                0.5*t(as.vector(exp(-xx_R1%*%Beta_1[,3]))*(Im(-y[1,2]*Conj(y[1,3]) + y[1,3]*Conj(y[1,2])))*t(xx_I1) +
                                              as.vector(exp(-xx_R1%*%Beta_1[,3]))*(2*abs(y[1,2])^2)*as.vector(xx_I1%*%Beta_2[,3])*t(xx_I1) +
                                              (1i)*as.vector(exp(-xx_R1%*%Beta_1[,3]))*(-y[1,1]*Conj(y[1,2])*(a[1]) + Conj(y[1,1]*Conj(y[1,2])*(a[1])))*t(xx_I1)) +
                                0.5*t(as.vector(exp(-xx_Re%*%Beta_1[,3]))*(Im(-y[eee,2]*Conj(y[eee,3]) + y[eee,3]*Conj(y[eee,2])))*t(xx_Ie) +
                                              as.vector(exp(-xx_Re%*%Beta_1[,3]))*(2*abs(y[eee,2])^2)*as.vector(xx_Ie%*%Beta_2[,3])*t(xx_Ie) +
                                              (1i)*as.vector(exp(-xx_Re%*%Beta_1[,3]))*(-y[eee,1]*Conj(y[eee,2])*(a[eee]) + Conj(y[eee,1]*Conj(y[eee,2])*(a[eee])))*t(xx_Ie))
                }
                gr <- c(gr1,gr2,gr3,gr4,gr5,gr6,gr7,gr8,gr9)
                gr_index <- (1:(dimen^2*nBeta))*kronecker(chol_index[Phi_temp,],rep(1,nBeta))
                gr_index <- gr_index[which(gr_index!=0)]
                gr <- Re(gr[gr_index])
                gr <- -gr
        }
        return(gr)
}