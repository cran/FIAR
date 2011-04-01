diffGranger <-
function(data,nx=1,ny=1,order=1,boot=FALSE,bs=100, p=.05){
data <- as.matrix(data)
X <- t(data)
X <- X-rowMeans(X)
data <- t(X)

if (boot == FALSE){
x <- as.matrix(data[,1:nx])

y <- as.matrix(data[,-(1:nx)][,1:ny])

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

Fxy = log(det(e1)/det(e2))

e1 = cov(x2)-cov(x2,cbind(xlag,zlag))%*%solve(cov(cbind(xlag,zlag)))%*%cov(cbind(xlag,zlag),x2)
e2 = cov(x2)-cov(x2,cbind(xlag,ylag,zlag))%*%solve(cov(cbind(xlag,ylag,zlag)))%*%cov(cbind(xlag,ylag,zlag),x2)

Fyx = log(det(e1)/det(e2))

	
      F <- Fxy - Fyx
      F
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

Fxy = log(det(e1)/det(e2))

e1 = cov(x2)-cov(x2,xlag)%*%solve(cov(xlag))%*%cov(xlag,x2)
e2 = cov(x2)-cov(x2,cbind(xlag,ylag))%*%solve(cov(cbind(xlag,ylag)))%*%cov(cbind(xlag,ylag),x2)

Fyx = log(det(e1)/det(e2))

	
      F <- Fxy - Fyx
      F
}
} 
 
 else{
 l <- median(b.star(data,round=TRUE)[,1]) # package {np}
 out <- tsboot(data,diffGranger,R=bs, l=l, sim ='geom',nx=nx,ny=ny,order=order) # package {boot}
 if(quantile(out$t-out$t0,1-p/ncol(data)/2)>out$t0){
out$sig=0} else {out$sig=1}
 out} 

}

