partGranger <-
function(data,nx=1,ny=1,order=1,boot=FALSE,bs=100,b=NA,p=.05){
data <- as.matrix(data)
X <- t(data)
X <- X-rowMeans(X)
data <- t(X)
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

e1 <- cov(y2)-cov(y2,cbind(ylag,zlag,z2))%*%solve(cov(cbind(ylag,zlag,z2)))%*%cov(cbind(ylag,zlag,z2),y2)
e2 <- cov(y2)-cov(y2,cbind(ylag,xlag,zlag,z2))%*%solve(cov(cbind(ylag,xlag,zlag,z2)))%*%cov(cbind(ylag,xlag,zlag,z2),y2)

F <- log(det(e1)/det(e2))
#fitF <-lm(y2~cbind(ylag,zlag,xlag,z2))
#fitR <- lm(y2~cbind(ylag,zlag,z2))
#SSEF <- sum(fitF$res^2)
#SSER <- sum(fitR$res^2)
#df1 <- order*nx
df2 <- nrow(x2)-order*(nx+ny+nz)-nz
dfr <- nrow(x2)-order*(ny+nz)-nz
#df2 <- (nrow(data)-order-(ncol(data)*order))
#  f <- ((SSER-SSEF)/df1)/(SSEF/df2)
prb <- 1-pf(exp(F),dfr,df2)
#prb <- 1-pf(f,df1,df2)
out <- list()
out$orig <- F
out$prob <- prb
out
}
if(boot==TRUE & (nx > 1 | ny > 1)){
out <-  " Warning: No bootstrap method available for Multivariate G-causality!! "
}
if(boot==TRUE & nx == 1 & ny == 1){
l <- ncol(data)
r <- nrow(data)
out<-list()
out$prob <- matrix(0,l,l)
out$orig <- matrix(0,l,l)
PG <- array(0,dim=c(l,l,bs))


for (bss in 1:bs){
XX  <- matrix(0,l,nrow(data)) 
nwin<- nrow(data)/b
inx <- rep(1:nwin,each=b)
ct  <-1
for (kk in 1:nwin){
temp<- sample(nwin)
winx<- which(inx==temp[1])
XX[,ct:(ct+b-1)]<-data[winx,]
ct <- ct + b
}
for(ii in 1:(l-1)) {
for (jj in (ii+1):l){
if (bss==1){
out$orig[ii,jj]=partGranger(cbind(data[,ii],data[,jj],data[,-c(ii,jj)]),nx=nx,ny=ny,order=order)$orig
out$orig[jj,ii]=partGranger(cbind(data[,jj],data[,ii],data[,-c(ii,jj)]),nx=nx,ny=ny,order=order)$orig
}

Xi <- rbind(XX[ii,],XX[jj,],XX[-c(ii,jj),])
Xir <- Xi[-((nx+1):(nx+ny)),] 
Exy<-armorf(t(Xi),nwin,b,order)
Ey<-armorf(t(Xir),nwin,b,order)

S11 <- Ey[1:ny,1:ny]
S22 <- Ey[(ny+1):ncol(Ey),(ny+1):ncol(Ey)]
S12 <- Ey[1:ny,(ny+1):ncol(Ey)]

Sig11 <- Exy[(1:ny),(1:ny)]
Sig22 <- Exy[((nx+ny+1):ncol(Exy)),(nx+ny+1):ncol(Exy)]
Sig12 <- Exy[(1:ny),(nx+ny+1):ncol(Exy)]

PG[ii,jj,bss] <- log(det(S11-S12%*%solve(S22)%*%S12) / det(Sig11-Sig12%*%solve(Sig22)%*%Sig12))-out$orig[ii,jj]

Xj <- rbind(XX[jj,],XX[ii,],XX[-c(ii,jj),])
Xjr <- Xj[-((nx+1):(nx+ny)),] 
Exy<-armorf(t(Xj),nwin,b,order)
Ey<-armorf(t(Xjr),nwin,b,order)

S11 <- Ey[1:ny,1:ny]
S22 <- Ey[(ny+1):ncol(Ey),(ny+1):ncol(Ey)]
S12 <- Ey[1:ny,(ny+1):ncol(Ey)]

Sig11 <- Exy[(1:ny),(1:ny)]
Sig22 <- Exy[((nx+ny+1):ncol(Exy)),(nx+ny+1):ncol(Exy)]
Sig12 <- Exy[(1:ny),(nx+ny+1):ncol(Exy)]

PG[jj,ii,bss] <- log(det(S11-S12%*%solve(S22)%*%S12) / det(Sig11-Sig12%*%solve(Sig22)%*%Sig12))-out$orig[jj,ii]

}}
}

for(ii in 1:(l-1)) {
for (jj in (ii+1):l){
if (quantile(PG[ii,jj,],1-p/ncol(data)/2)>out$orig[ii,jj])
{out$prob[ii,jj]=0} else {out$prob[ii,jj]=1}

if (quantile(PG[jj,ii,],1-p/ncol(data)/2)>out$orig[jj,ii])
{out$prob[jj,ii]=0} else {out$prob[jj,ii]=1}

}}
}
out
}
