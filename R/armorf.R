armorf <-  function(X,Nr,Nl,p){


X <- t(as.matrix(X))
N <- ncol(X)
L <- nrow(X)
X <- X-rowMeans(X)

 
  startinx <- 1
  for (i in 1:Nr){
  endinx <- startinx + Nl - 1
   if (endinx > N) endinx <- N
   XX <- X[,startinx:endinx]
   XX <- XX - rowMeans(XX)
   X[,startinx:endinx] <- XX
   startinx <- startinx + Nl   
                           }
x <- X                           
R0 <- matrix(0,L,L)
R0f <- R0
R0b <- R0
pf <- R0
pb <- R0
pfb <- R0
En <- R0
ap <- array(0,dim <- c(L,L,p+1))
bp <- ap
a  <-  array(0,dim <- c(L,L,p+1))
b  <-  a

ap[,,1] <- R0
bp[,,1] <- R0
for (i in 1:Nr){
En <- En + x[,((i-1)*Nl+1):(i*Nl)]%*%t(x[,((i-1)*Nl+1):(i*Nl)]) 
ap[,,1] <- ap[,,1] + x[,((i-1)*Nl+2):(i*Nl)] %*% t(x[,((i-1)*Nl+2):(i*Nl)])
bp[,,1] <- bp[,,1] + x[,((i-1)*Nl+1):(i*Nl-1)] %*% t(x[,((i-1)*Nl+1):(i*Nl-1)])
}

ap[,,1] <- solve(t(chol(ap[,,1]/Nr*(Nl-1))))
bp[,,1] <- solve(t(chol(bp[,,1]/Nr*(Nl-1))))

for (i in 1:Nr){
efp <- ap[,,1]%*% x[,((i-1)*Nl+2):(i*Nl)]
ebp <- bp[,,1]%*% x[,((i-1)*Nl+1):(i*Nl-1)]
pf  <-  pf + efp%*%t(efp)
pb  <-  pb + ebp%*%t(ebp)
pfb  <-  pfb + efp%*%t(ebp)
}

En  <-  t(chol(En/N))
## initialize
#coeff  <-  logical(0)
kr  <-  logical(0)

for (m in 1:p){
## Calculate the next order reflection (parcor) coefficient
ck  <-  solve(t(chol(pf)))%*%pfb%*%solve(chol(pb))
kr  <-  cbind(kr,ck)
## Update the forward and backward prediction errors
ef  <-  diag(L)-ck%*%t(ck)
eb  <-  diag(L)-t(ck)%*%ck

## Update the prediction error

En  <-  En%*%t(chol(ef))
E  <-  (ef+eb)/2

## Update the coefficients of the forward and backward prediction errors
#ap[,,m+1] <- diag(L)-diag(L)
#bp[,,m+1] <- diag(L)-diag(L)
pf  <-  diag(L)-diag(L)
pb  <-  diag(L)-diag(L)
pfb  <-  diag(L)-diag(L)

for (i in 1:(m+1)){
a[,,i]  <-  solve(t(chol(ef)))%*%(ap[,,i]-ck%*%bp[,,(m+2-i)])
b[,,i]  <-  solve(t(chol(eb)))%*%(bp[,,i]-t(ck)%*%ap[,,(m+2-i)])
}
for(k in 1:Nr){
efp  <-  matrix(0,L,(Nl-m-1))
ebp  <-  matrix(0,L,(Nl-m-1))

  for (i in 1:(m+1)){
k1  <-  m+2-i+(k-1)*Nl+1
k2  <-  Nl-i+1+(k-1)*Nl
efp  <-  efp + a[,,i]%*%x[,k1:k2]
ebp  <-  ebp + b[,,m+2-i]%*%x[,(k1-1):(k2-1)]
}
pf  <-  pf + efp%*%t(efp)
pb  <-  pb + ebp%*%t(ebp)
pfb <-  pfb + efp%*%t(ebp)
}
ap[,,1:(m+1)]  <-  a[,,1:(m+1)]
bp[,,1:(m+1)]  <-  b[,,1:(m+1)]

}
E  <-  En %*% t(En)
E
}

