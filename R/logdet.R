
logdet<-function(C){
# if(abs(log(det(as.matrix(S))))==Inf){
# sv <- sqrt(abs(spm_eigen(as.array(tcrossprod(S)))))
# #sv<-sqrt(sv)
# Sv <- sum(log(sv[sv>10^-16 && sv<(1/10^-16)]))
# Sv
#                                     }else{
# Sv=log(det(as.matrix(S)))
# Sv
#      }
# }

TOL=1e-16
n = dim(C)[1]
s=diag(C)
i=which(s>TOL & s<(1/TOL))
C=C[i,i]
H=sum(log(diag(C)))
if (abs(H)==Inf){
s  = sqrt(abs(eigen(as.array(tcrossprod(C)))$values))
    H  = sum(log(s(s > TOL & s < 1/TOL)))
H} else{H=H}
H}
