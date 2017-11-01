Beta_derive2 <- function(x, yobs_tmp, chol_index, Phi_temp,
                         tau_temp_1, tau_temp_2, Beta_temp_1, Beta_temp_2,
                         sigmasqalpha, nbasis, nseg){
        
        eee <- dim(yobs_tmp)[1]
        yobs_tmp_1 <- yobs_tmp[1:nseg,]
        yobs_tmp_2 <- yobs_tmp[(nseg+1):eee,]
        output1 <- Beta_derive1(x, yobs_tmp_1, chol_index, Phi_temp, 
                                tau_temp_1, Beta_temp_1, sigmasqalpha, nbasis, global_info)
        f1 <- output1$value
        grad1 <- output1$gradient
        hes1 <- output1$hessian
        
        output2 <- Beta_derive1(x, yobs_tmp_2, chol_index, Phi_temp, 
                                tau_temp_2, Beta_temp_2, sigmasqalpha, nbasis, global_info)
        f2 <- output2$value
        grad2 <- output2$gradient
        hes2 <- output2$hessian

        f <- f1 + f2
        gr <- grad1 + grad2
        h <- hes1 + hes2
        list(value = f, gradient = gr, hessian = h)
}