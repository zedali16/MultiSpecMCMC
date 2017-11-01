Gradient2 <- function(yobs_tmp, chol_index, Phi_temp, tau_temp_1, tau_temp_2,
                         Beta_temp_1, Beta_temp_2, sigmasqalpha, nbasis, nseg, global_info){
        
        dimen <- global_info$dimen
        yobs_tmp_1 <- yobs_tmp[1:nseg,]
        yobs_tmp_2 <- yobs_tmp[-(1:nseg),]
        grad1 <- Gradient1(yobs_tmp_1, chol_index, Phi_temp, tau_temp_1, Beta_temp_1, sigmasqalpha, nbasis, global_info)
        grad2 <- Gradient1(yobs_tmp_2, chol_index, Phi_temp, tau_temp_2, Beta_temp_2, sigmasqalpha, nbasis, global_info)
        gr <- grad1 + grad2
        return(gr)
}