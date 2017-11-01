nchoosek <- function (n, k){
                # function for calculating all subsets of size k from n objects
                # taken from package vsn (provided by Wolfgang Huber under LGPL)
                # slightly modified to also work for n=k=2 by Ulrike Groemping
                if (!is.numeric(n) || !is.numeric(k) || is.na(n) || is.na(k) || 
                    length(n) != 1 || length(k) != 1) 
                        stop("arguments must be non-NA numeric scalars.")
                if (k > n || k < 0) 
                        stop("Arguments must satisfy 0 <= k <= n.")
                nck = choose(n, k)
                res = matrix(NA, nrow = k, ncol = nck)
                res[, 1] = 1:k
                j = 2
                repeat {
                        if (j > nck) 
                                break
                        res[, j] = res[, j - 1]
                        i = k
                        repeat {
                                res[i, j] = res[i, j] + 1
                                if (res[i, j] <= n - (k - i)) 
                                        break
                                i = i - 1
                                stopifnot(i >= 1)
                        }
                        if (i < k) 
                                res[(i + 1):k, j] = res[i, j] + 1:(k - i)
                        j = j + 1
                }
                stopifnot(all(res[, nck] == (n - k + 1):n))
                stopifnot(all(res <= n) && all(res >= 1))
                return(res)
}


chol_ind <- function(dimen){
        chol_index = matrix(0,2^(dimen^2),dimen^2)
        k=0
        for(i in 0:dimen^2){
                C = nchoosek(dimen^2,i)
                if(i==0){
                        k = k+1;
                        chol_index[1,] = 0
                } else{
                        for(j in 1:dim(C)[2]){
                                k = k+1
                                chol_index[k,C[,j]] = 1
                        }
                }
        }
        return(chol_index)
}

##################################################
#                  matrix power                  #
##################################################
matpower <- function(a,alpha){
        small <- .000000001
        if (length(c(a))==1){
                ai=a^alpha
        }else{
                p1<-nrow(a)
                eva<-eigen(a)$values
                eve<-eigen(a)$vectors
                eve<-eve/t(matrix((diag(t(eve)%*%eve)^0.5),p1,p1))
                index<-(1:p1)[eva>small]
                evai<-eva
                evai[index]<-(eva[index])^(alpha)
                ai<-eve%*%diag(evai)%*%t(eve)
        }
        return(ai)
}