lin_basis_func <- function(freq_hat,global_info){
        nBeta = global_info$nBeta
        nfreq_hat = length(freq_hat)
        xx_r = matrix(1,nfreq_hat,nBeta)
        xx_i = matrix(1,nfreq_hat,nBeta)
        for(j in 2:nBeta){
                xx_r[,j] = sqrt(2)*cos(2*pi*(j-1)*freq_hat)/(2*pi*(j-1))
        }
        for(j in 1:nBeta){
                xx_i[,j] = sqrt(2)*sin(2*pi*(j)*freq_hat)/(2*pi*(j))
        }
        return(list(xx_r=xx_r, xx_i=xx_i))
}