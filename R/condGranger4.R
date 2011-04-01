condGranger2 <-
function(data,nx=1,ny=1,order=1,boot=FALSE,prob=TRUE,bs=100, p=.05){

l <- ncol(data)
r <- nrow(data)
out<-list()
out$prb <- matrix(0,l,l)
out$orig <- matrix(0,l,l)
cg <- array(0,dim=c(l,l,bs))
ll <- median(b.star(data,round=TRUE)[,1])

for (bss in 1:bs){
XX  <- matrix(0,l,r) 
nwin<- floor(r/ll)
inx <- rep(1:nwin,each=ll)
ct  <-1
for (kk in 1:nwin){
temp<- sample(nwin)
winx<- which(inx==temp[1])
XX[,ct:(ct+ll-1)]<-data[winx,]
ct <- ct + ll
}
XX <- XX[,1:max(which(apply(XX,2,prod)!=0))]
for(ii in 1:(l-1)) {
for (jj in (ii+1):l){
if (bss==1){
out$orig[ii,jj]=condGranger(cbind(data[,ii],data[,jj],data[,-c(ii,jj)]),nx=nx,ny=ny,order=order)$orig
out$orig[jj,ii]=condGranger(cbind(data[,jj],data[,ii],data[,-c(ii,jj)]),nx=nx,ny=ny,order=order)$orig
}

Xi <- t(rbind(XX[ii,],XX[jj,],XX[-c(ii,jj),]))
cg[ii,jj,bss]=  condGranger(Xi,nx=1,ny=1,order=order)$orig
Xj <- t(rbind(XX[jj,],XX[ii,],XX[-c(ii,jj),]))
cg[jj,ii,bss]=  condGranger(Xj,nx=1,ny=1,order=order)$orig
}}
}
for(ii in 1:(l-1)) {
for (jj in (ii+1):l){
if (quantile(cg[ii,jj,],1-p/(ncol(data)^2-ncol(data)))>out$orig[ii,jj])
{out$prb[jj,ii]=0} else {out$prb[jj,ii]=1}

if (quantile(cg[jj,ii,],1-p/(ncol(data)^2-ncol(data)))>out$orig[jj,ii])
{out$prb[ii,jj]=0} else {out$prb[ii,jj]=1}

}}
out$orig<-t(out$orig)
out
}
