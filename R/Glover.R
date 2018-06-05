
Glover <- function(DCM,ts){

Kernel <- spm_kernels(DCM)

na <- matrix(0,DCM$v,DCM$n)

for (i in 1:DCM$n){
H <- fft(Kernel[,i,1])
M <- fft(ts[,i])
N0 <- DCM$v*DCM$Ce[(i-1)*DCM$v+1,(i-1)*DCM$v+1]
na[,i] <- as.double(as.double(fft(Conj(H)*M)/(as.double(H*Conj(H))+N0),inverse=TRUE))/DCM$v
#na[,i] <- as.real(as.real(fft(Conj(H)*M)/(as.real(H*Conj(H))+N0),inverse=TRUE))/DCM$v
#na[,i] <- as.real(fft(Conj(H)*M/(as.real(H*Conj(H))),inverse=TRUE))/DCM$v
}
#na<-as.matrix(na)
na
}
