
spm_int1 <-
function(DCM){
       x_exp <- rbind(1,as.matrix(rep(0,DCM$x)))
       #PM=c(0,DCM$a,DCM$b,DCM$c,DCM$H%x%rep(1,DCM$n))

    # Integrate
      bireduce<-spm_bireduce(DCM$Ep,DCM$x,DCM$n,DCM$m,DCM$TE)
        M0 <- bireduce$M0
        M1 <- bireduce$M1

   # Generate Bold signal
    y <- matrix(0,DCM$n,DCM$v)
    dy <- y
    ns <- logical()
      for (i in 1:length(DCM$T0)){
        u <- DCM$sf[DCM$T0[i],]
        if (DCM$s[i]>DCM$v){

        J <- M0
        for (j in 1:DCM$m){
        J <- J+u[j]*M1[,,j]
                      }
                       }
          else{y[,DCM$s[i]] <- funl2(x_exp[2:length(x_exp)], DCM$Ep,DCM$n,DCM$TE)

                  ns <- cbind(ns,x_exp)
}


      x_exp <- spm_expm2((J*DCM$dt0[i]),x_exp)

      }

  ts <-  t(ns[2:(DCM$n+1),])
  #ts <-  t(ns[(1*DCM$n+2):(2*DCM$n+1),])
  #ts <-  t(ns[(2*DCM$n+2):(3*DCM$n+1),])
  #ts <-  t(ns[(3*DCM$n+2):(4*DCM$n+1),])
  #ts <-  t(ns[(4*DCM$n+2):(5*DCM$n+1),])
 }
