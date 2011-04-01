condGranger <-
function(data,nx=1,ny=1,order=1,boot=FALSE,prob=TRUE,bs=100, p=.05){
data <- as.matrix(data)
X <- t(data)
X <- X-rowMeans(X)
data <- t(X)

  if (boot == FALSE){
     x <- as.matrix(data[,1:nx])

     y <- as.matrix(data[,-(1:nx)])[,1:ny]

       if (sum(as.matrix(data[,-(1:(nx+ny))]))!=0)
         {
           z <- as.matrix(data[,-(1:(nx+ny))])
           nz <- ncol(z)

# Compute F2

                x  <-  embed(x,order+1)
                x2 <- as.matrix(x[,1:nx])
                xlag <- as.matrix(x[,-(1:nx)])

                y  <-  embed(y,order+1)
                y2 <- as.matrix(y[,1:ny])
                ylag <- as.matrix(y[,-(1:ny)])

                z  <-  embed(z,order+1)
                z2 <- as.matrix(z[,1:nz])
                zlag <- as.matrix(z[,-(1:nz)])

             e1 = cov(y2)-cov(y2,cbind(ylag,zlag))%*%solve(cov(cbind(ylag,zlag)))%*%cov(cbind(ylag,zlag),y2)
             e2 = cov(y2)-cov(y2,cbind(ylag,xlag,zlag))%*%solve(cov(cbind(ylag,xlag,zlag)))%*%cov(cbind(ylag,xlag,zlag),y2)

           F = log(det(e1)/det(e2))

	    fitx <-lm(y2~0+cbind(ylag,zlag,xlag))
	    fit <- lm(y2~0+cbind(ylag,zlag))
	    RSS0 <- sum(fit$res^2)
	    RSS1 <- sum(fitx$res^2)
	    df1 <- order*nx
	    df2 <- nrow(x2)-order*(nx+ny+nz)
	    dfr <- nrow(x2)-order*(ny+nz)
	    f <- ((RSS0-RSS1)/df1)/(RSS1/df2)
	    #prb <- 1-pf(exp(F),order,df2)
	    prb <- 1-pf(exp(F),dfr,df2)
	    if(prob==TRUE){out <- list()
	    out$orig <- F
	    out$prob <- prb}
	    else{out <- F}
	    out
	    }
else{
x  <-  embed(x,order+1)
x2 <- as.matrix(x[,1:nx])
xlag <- as.matrix(x[,-(1:nx)])

y  <-  embed(y,order+1)
y2 <- as.matrix(y[,1:ny])
ylag <- as.matrix(y[,-(1:ny)])

e1 = cov(y2)-cov(y2,ylag)%*%solve(cov(ylag))%*%cov(ylag,y2)
e2 = cov(y2)-cov(y2,cbind(ylag,xlag))%*%solve(cov(cbind(ylag,xlag)))%*%cov(cbind(ylag,xlag),y2)

F = log(det(e1)/det(e2))

fitx <-lm(y2~cbind(ylag,xlag))
fit <- lm(y2~cbind(ylag))
RSS0 <- sum(fit$res^2)
RSS1 <- sum(fitx$res^2)
df1 <- order*nx
df2 <- nrow(x2)-order*(ny+nx)
dfr <- nrow(x2)-order*(ny)
f <- ((RSS0-RSS1)/df1)/(RSS1/df2)
#prb <- 1-pf(exp(F),dfr,df2)
if(prob==TRUE){out <- list()
out$orig <- f
out$prob <- prb}
else{out <- f}
out
}
}
else{
 l <- median(b.star(data,round=TRUE)[,1]) # package {np}
 out <- tsboot(data,condGranger,R=bs,l=l, sim ='geom',nx=nx,ny=ny,order=order,prob=FALSE) # package {boot}
 if(quantile(out$t-out$t0,1-p/ncol(data)/2)>out$t0){
 out$sig=0} else {out$sig=1}
 out
}
}
