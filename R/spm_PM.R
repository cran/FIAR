spm_PM <-
function(a,b,c,h,n,m){
	
	a <- matrix(a,ncol=n,byrow=TRUE)
        diag(a) <- -1
	dim(b) <- c(n,n,m)
                         for (i in 1:m){
                              b[,,i] <- t(b[,,i])
                                        }
		PM=c(0,a,b,c,h%x%rep(1,n))
		PM
		                 }

