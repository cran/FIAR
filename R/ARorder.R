ARorder <- function (data, min = 1, max = 10,type='AIC') 
{
if (type=='AIC'){
   
 AIC <- rep(Inf, max)
    d <- as.matrix(data)
    N <- NCOL(d)
        if(N > 1){
	      for (i in min:max) {
		z <- embed(d, i + 1)
		z2 <- z[, (1:N)]
		zlag <- z[, -(1:N)]
		fit <- lm(z2 ~ zlag)
		AIC[i] <- log(det(cov(fit$res))) + 2 * N^2 * i/NROW(d)
	                          }
    orr <- which(AIC == min(AIC[abs(AIC) != Inf]))
    
                  }
        else{
	      for (i in min:max) {
		z <- embed(d, i + 1)
		z2 <- z[, (1:N)]
		zlag <- z[, -(1:N)]
		fit <- lm(z2 ~ zlag)
		AIC[i] <- log(var(fit$res)) + 2 * N^2 * i/NROW(d)
	      }
             orr <- which(AIC == min(AIC[abs(AIC) != Inf]))
    
             }
} 
if (type=='BIC'){
BIC <- rep(Inf,max)
    d <- as.matrix(data)
    N <- NCOL(d)
        if(N > 1){
           for (i in min:max){
             z  <-  embed(d,i+1)
             z2 <- z[,(1:N)]
             zlag <- z[,-(1:N)]
             fit <- lm(z2~zlag)
             BIC[i] <- log(det(cov(fit$res)))+log(NROW(d))*N^2*i/NROW(d)
                             }
          orr <- which(BIC==min(BIC[abs(BIC)!=Inf]))
    
                 }

	else{
	    for (i in min:max) {
	      z <- embed(d, i + 1)
	      z2 <- z[, (1:N)]
	      zlag <- z[, -(1:N)]
	      fit <- lm(z2 ~ zlag)
	      BIC[i] <- log(var(fit$res)) + log(NROW(d))*N^2*i/NROW(d)
	    }
	  orr <- which(BIC == min(BIC[abs(BIC) != Inf]))

             }

}
orr
}
