MultiSpect_partition <- function(zt,output){
        attach(output)
        nobs = dim(zt)[1]; dimen=dim(zt)[2]
        dev.new()
        hist(nexp_curr[(nwarmup+1):nloop],breaks=0:nexp_max,xlab="Partitions",main="Histogram of Number of Partitions")
        for (j in 1:nexp_max){
                kk <- which(nexp_curr[(nwarmup + 1):nloop] == j)
                if (length(kk) != 0 & j > 1){
                        for (k in 1:(j - 1)){
                                dev.new()
                                plot(xi[[j]][k, kk + nwarmup], type = "l", ylab="Time Point", xlab="Index",
                                     main = paste("Plot of Location of Partition Point", k, 'Given', j, 'Segments'))
                        }
                        for (k in 1:(j - 1)){
                                dev.new()
                                hist(xi[[j]][k, kk + nwarmup],1:nobs,xlab="Time Point",
                                     main=paste("Histogram of Location of Partition Point", k, 'Given', j, 'Segments'))
                        }
                }
        }
        
        prob <- histc(nexp_curr[(nwarmup+1):nloop],1:nexp_max)$cnt/length((nwarmup+1):nloop)
        detach(output)
        return(prob)
}
