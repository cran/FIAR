
spm_expm2 <-
function(J,x){
    x0 <- x
    fx <- (J)%*%x
    j  <- 1

	#while (eigen2(tcrossprod(fx))$values[1] > 1e-16){
#while (spm_eigen(tcrossprod(fx))[NROW(fx)] > 1e-16){
    while (eigen(tcrossprod(fx))$values[1] > 1e-16){
 	j  <- j + 1
	x  <- x + fx
	fx <- (J)%*%fx/j

  	if

	#(eigen2(tcrossprod(x))$values[1] > 1e16){
	#(spm_eigen(tcrossprod(x))[NROW(x)] > 1e16){
	(eigen(tcrossprod(x))$values[1] > 1e16){
	x <- spm_expm1((J))%*%x0}
	}
    x
    }

