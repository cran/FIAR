#' @export
dcmGenerate <-
function(DCM=DCM,SNR=0,ar=0,names=DCM$names){
    PM <- spm_PM(DCM$a,DCM$b,DCM$c,DCM$h,DCM$n,DCM$m)
    #X0 <- HPF(DCM)
    DCM<-stimfun(DCM)
    y <- spm_int(PM,DCM)
    y <- t(y)
      # add noise
    if(SNR>0){r=diag(apply(t(y),2,sd)/SNR)
    p       = 1
    #a       = 0    # AR(1) coeff: for the moment set to zero
    a <- diag(DCM$v)
    for (i in 2:DCM$v){
    a[i,(i-1)]<--ar
    }
    K=solve(a)
    K=K*sqrt(DCM$v/sum(diag(tcrossprod(K))))
    z=matrix(rnorm(DCM$n*DCM$v),DCM$v,DCM$n)

    e=K%*%z
    y=y+t(e%*%r)
	    }

    y <- t(y)

#	Now orthogonalise data with respect to effects of no interest
#   If X0 is just a vector of 1s this amounts to making the data zero mean

	X0 <- rep(1,DCM$v)

	Xp  <- tcrossprod((X0%*%solve(crossprod(X0))),X0)

	for (i in 1:DCM$n){
    y[,i] <- y[,i]-Xp%*%y[,i]
	                  }
if(length(names)==ncol(y)) colnames(y)<-names
y
	}

