pdiffGranger <-
function(data,nx=1,ny=1,order=1,boot=FALSE,bs=100){
data <- as.matrix(data)
if (boot == FALSE){
x <- as.matrix(data[,1:nx])

y <- as.matrix(data[,-(1:nx)][,1:ny])


z <- as.matrix(data[,-(1:(nx+ny))])
nz <- ncol(z)

# Compute F1

x  <-  embed(x,order+1)
x2 <- as.matrix(x[,1:nx])
xlag <- as.matrix(x[,-(1:nx)])

y  <-  embed(y,order+1)
y2 <- as.matrix(y[,1:ny])
ylag <- as.matrix(y[,-(1:ny)])

z  <-  embed(z,order+1)
z2 <- as.matrix(z[,1:nz])
zlag <- as.matrix(z[,-(1:nz)])

e1 = cov(y2)-cov(y2,cbind(ylag,zlag,z2))%*%solve(cov(cbind(ylag,zlag,z2)))%*%cov(cbind(ylag,zlag,z2),y2)
e2 = cov(y2)-cov(y2,cbind(ylag,xlag,zlag,z2))%*%solve(cov(cbind(ylag,xlag,zlag,z2)))%*%cov(cbind(ylag,xlag,zlag,z2),y2)

Fxy = log(det(e1)/det(e2))

e1 = cov(x2)-cov(x2,cbind(xlag,zlag,z2))%*%solve(cov(cbind(xlag,zlag,z2)))%*%cov(cbind(xlag,zlag,z2),x2)
e2 = cov(x2)-cov(x2,cbind(xlag,ylag,zlag,z2))%*%solve(cov(cbind(xlag,ylag,zlag,z2)))%*%cov(cbind(xlag,ylag,zlag,z2),x2)

Fyx = log(det(e1)/det(e2))

	
      F <- Fxy - Fyx
      F
 }
 else{
 l <- median(b.star(data,round=TRUE)[,1]) # package {np}
 out <- tsboot(data,pdiffGranger,R=bs, l=l, sim ='fixed',nx=nx,ny=ny,order=order) # package {boot}
 out  }
}

