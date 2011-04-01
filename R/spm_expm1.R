spm_expm1 <-
function(J){
    e= frexp(norm(J,'I'))[2]
    s     = max(0,e+1);
    J <- (J)/2^s
    X <- J
    c <- .5
    E <- diag(sqrt(length(J)))+c*J
    D <- diag(sqrt(length(J)))-c*J
    q <- 6
    p <- 1
	for (k in 2:q){
	c   <- c * (q-k+1) / (k*(2*q-k+1))
	X   <- J%*%X
	cX  <- c*X
	E   <- E + cX
	if (p==1)(D =D + cX)
	if (p!=1)(D = D - cX)
	p =!p
	}
    E <- solve(D)%*%E
    for (k in 1:s){
    E <- E%*%E
    }
    E
    }

