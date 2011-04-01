funx <-
function(x,u0=0,PM,n,m){
    q <- exp(PM[1])
PM <- PM[-1]
j <- 1:(n*n)
A0 <- matrix(PM[j],n,n)*q
PM <- PM[-j]
j <- 1:(n*n*m)
B0 <- PM[j]
dim(B0) <- c(n,n,m)
for (i in 1:m){
	B0[,,i] <- B0[,,i]*q
              }
PM <- PM[-j]
j <- 1:(n*m)
C0 <- matrix(PM[j],n,m)
PM <- PM[-j]
H0 <- matrix(PM,nrow=n)

	
	for(i in 1:m){
	A0 <- A0+u0[i]*B0[,,i]}
    x <- matrix(x,nrow=n)
    x[,2:5] <- exp(x[,2:5])
    x[,2] <- x[,2]-1
    fv <- x[,4]^(1/H0[,4])
    ff <- (1 - (1 - H0[,5])^(1/x[,3]))/H0[,5];
    y         <- matrix(0,n,5);
    y[,1]    <- A0%*%x[,1] + C0%*%u0;
    y[,2]    <- (x[,1]  - H0[,1]*x[,2] - H0[,2]*(x[,3] - 1))/(x[,2] + 1);
    y[,3]    <- x[,2]/x[,3];
    y[,4]    <- (x[,3] - fv)/(H0[,3]*x[,4]);
    y[,5]    <- (ff*x[,3] - fv*x[,5]/x[,4])/(H0[,3]*x[,5]);
    y <- matrix(y,ncol=1)
    y
    }


#ts1000=dcmGenerate(DCM)
