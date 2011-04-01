spm_kernels <- function(DCM){

bireduce<-spm_bireduce(DCM$Ep,DCM$x,DCM$n,DCM$m,DCM$TE)

dt <- DCM$TR/DCM$T
M0 <- bireduce$M0
M1 <- bireduce$M1
L1 <- bireduce$L1
n <- ncol(M0)
m <- DCM$m
l <- DCM$n
N <- DCM$v
H1 <- array(0,c(N,n,m))
K1 <- array(0,c(N,l,m))
K2 <- array(0,c(N,N,l,m,m))

e1 <- expm(dt*M0)
e2 <- expm(-dt*M0)

M <- array(0,c(dim(M1),N))

for (p in 1:m){

      M[,,p,1] <- as.matrix(e1%*%M1[,,p]%*%e2)
              }

for (i in 2:N){
     for (p in 1:m){
           M[,,p,i] <- as.matrix(e1%*%M[,,p,i-1]%*%e2)
                   }
              }

X0 <- as.vector(c(1,rep(0,n-1)))

 mp <- function(mat,pow){
 ans <- mat
 for ( i in 1:(pow-1)){
 ans <- mat%*%ans
 }
 return(ans)
 }

e1N <- mp(e1,N)
H0 <- e1N%*%X0
K0 <- L1%*%H0
H1 <- array(0,c(N,n,p))
K1 <- array(0,c(N,l,p))
for (p in 1:m){

   for (i in 1:N){

      H1[i,,p] <- as.matrix(M[,,p,i]%*%H0)
      K1[i,,p] <- as.matrix(H1[i,,p]%*%t(L1))

                 }
              }
K1
}
