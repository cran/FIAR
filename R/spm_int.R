
spm_int <-
function(PM,DCM){
       x_exp <- rbind(1,as.matrix(rep(0,DCM$x)))
       #PM=c(0,DCM$a,DCM$b,DCM$c,DCM$H%x%rep(1,DCM$n))

    # Integrate
      bireduce<-spm_bireduce(PM,DCM$x,DCM$n,DCM$m,DCM$TE)
        M0 <- bireduce$M0
        M1 <- bireduce$M1

   # Generate Bold signal
    y <- matrix(0,DCM$n,DCM$v)
    dy <- y

      for (i in 1:length(DCM$T0)){
	u <- DCM$sf[DCM$T0[i],]
	if (DCM$s[i]>DCM$v){

	J <- M0
	for (j in 1:DCM$m){
	J <- J+u[j]*M1[,,j]
                      }
	               }
	  else{y[,DCM$s[i]] <- funl2(x_exp[2:length(x_exp)], PM,DCM$n,DCM$TE)
                  #ns <- logical()
                  #ns <- cbind(ns,x_exp)
}


      x_exp <- spm_expm2((J*DCM$dt0[i]),x_exp)

      }

   t(y)
 }

