#' @export
diffGranger <-
function(data,nx=1,ny=1,order=1,perm=FALSE,bs=100){
data <- as.matrix(data)
X <- t(data)
X <- X-rowMeans(X)
data <- t(X)

  if (perm == FALSE){
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


	    fitF <-lm(y2~cbind(ylag,zlag,xlag))
	    fitR <- lm(y2~cbind(ylag,zlag))
	    SSER <- sum(fitR$res^2)
	    SSEF <- sum(fitF$res^2)
            fxy <- log(SSER/SSEF)
            fitF <-lm(x2~cbind(ylag,zlag,xlag))
	    fitR <- lm(x2~cbind(xlag,zlag))
	    SSER <- sum(fitR$res^2)
	    SSEF <- sum(fitF$res^2)
            fyx <- log(SSER/SSEF)
	    Fdiff <- fxy-fyx
            Fdiff
	    }
		  else{
		  x  <-  embed(x,order+1)
		  x2 <- as.matrix(x[,1:nx])
		  xlag <- as.matrix(x[,-(1:nx)])

		  y  <-  embed(y,order+1)
		  y2 <- as.matrix(y[,1:ny])
		  ylag <- as.matrix(y[,-(1:ny)])

		  fitF <-lm(y2~cbind(ylag,xlag))
		  fitR <- lm(y2~cbind(ylag))
		  SSER <- sum(fitR$res^2)
		  SSEF <- sum(fitF$res^2)
		  fxy <- log(SSER/SSEF)
		  fitF <-lm(x2~cbind(ylag,xlag))
		  fitR <- lm(x2~cbind(xlag))
		  SSER <- sum(fitR$res^2)
		  SSEF <- sum(fitF$res^2)
		  fyx <- log(SSER/SSEF)
                  Fdiff <- fxy-fyx
                  Fdiff
		  }
                          }


else{
l <- ncol(data)
r <- nrow(data)
out<-list()
cg <- matrix(0,bs,1)
ll <- b.star(data,round=TRUE)[,1]


for (bss in 1:bs){
XX  <- matrix(0,l,r)
  for (pp in 1:l){

    nwin<- floor(r/ll[pp])
    temp<- sample(nwin)
    inx <- rep(temp,each=ll[pp])
    data_ind <- cbind(inx,data[1:(ll[pp]*nwin),pp])
    XX[pp,1:(ll[pp]*nwin)] <- data_ind[order(data_ind[,1]),-1]

                }
  XX <- XX[,1:max(which(apply(XX,2,prod)!=0))]

if (bss==1){
out$orig  <- diffGranger(data,nx=nx,ny=ny,order=order)
           }
cg[bss,]    <- diffGranger(t(XX),nx=nx,ny=ny,order=order)

}
#if (quantile(cg,1-p/(ncol(data)^2-ncol(data)))>out$orig)
#{out$prb=0} else {out$prb=1}
if(out$orig>0){
out$prob <- length(cg[cg>=out$orig])/length(cg)}
if(out$orig<0){
out$prob <- length(cg[cg<=out$orig])/length(cg)}
out
}
}
