spm_dctmtx <-
function(N, K) {

    n <- 0:(N-1)    
    C <- matrix(0, nrow=N, ncol=K)

    C[,1] <- 1/sqrt(N)
    for(k in 2:K) {
        C[,k] = sqrt(2/N)*cos(pi*(2*n+1)*(k-1)/(2*N))
    }
 
    C
    }

