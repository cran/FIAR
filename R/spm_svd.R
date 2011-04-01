spm_svd <-
function(X,tol){
    X <- as.matrix(X)
    M <- dim(X)[1]
    N <- dim(X)[2]
    p <- which(rowSums(X)!=0)
    q <- which(colSums(X)!=0)
    X <- X[p,q]
    i <- which(X!=0,arr.ind=TRUE)[,1]
    j <- which(X!=0,arr.ind=TRUE)[,2]
    
    s <- X[X!=0]
    M_ <- dim(X)[1]
    N_ <- dim(X)[2]
    ve <- svd(X)$v
    uu <- svd(X)$u
    S <- diag(svd(X)$d)
   s <- diag(S)^2
    j <- which(s*length(s)/sum(s) >= tol & s >= 0)
    ve <- ve[,j]
    uu <- uu[,j]
    S <- S[j,j]
    j <- length(j)
    U <- matrix(0,M,j)
    V <- matrix(0,N,j)
    
    V[q,] <- ve
    V
    }

