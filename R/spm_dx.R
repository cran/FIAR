spm_dx <-
function(dfdx,f,t){
dfdx<-as.matrix(dfdx)
f <- as.matrix(f)
    #t<-t/sqrt(eigen2(dfdx%*%t(dfdx))$values[1])
    t<-t/sqrt(eigen(dfdx%*%t(dfdx))$values[1])
    if (t > 10^8){
	    dx = -solve(dfdx)%*%f
      }
    Jx=cbind(rbind(0,f),rbind(0,dfdx))
    dx=spm_expm1((Jx*t))
    dx=dx[2:length(dx[1,]),1]
    dx
    }

