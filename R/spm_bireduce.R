spm_bireduce <-
function(PM,x,n,m,TE){
    
    X <- matrix(0,x,1)
    u0 <- rep(0,m)
    Vx <- diag(x)
    Vu <- diag(m)
    f0 <- funx(X,u0=u0, PM,n,m)
    l0 <- funl2(X, PM,n,TE)
      # differentiate wrt x (dfdx) 
	dx <- exp(-8)
	dfdx <- matrix(0,x,x)
	
	    for (i in 1:x)
	    {
	    xmi <- X+Vx[,i]*dx
	    fi <- funx(xmi,u0=u0, PM,n,m)
	    dfdx[,i]  <- spm_dfdx(fi,f0,dx)
	    }
      
      # differentiate wrt x and u (dfdxu)
	
	
	dfdxu <- rep(NA,length(dfdx)*m)
    dim(dfdxu) <- c(dim(dfdx),m)
    
      for (j in 1:m){
    
         xmu <- u0 + Vu[,j]*dx
       J <- matrix(0,x,x)
    f0 <- funx(X,u0=xmu, PM,n,m)
    
        for (i in 1:x)
            {
        xmi <- X+Vx[,i]*dx
        fi <- funx(xmi,u0=xmu, PM,n,m)
        J[,i]  <- spm_dfdx(fi,f0,dx)
            }
        dfdxu[,,j] <- spm_dfdx(J,dfdx,dx)
                        } 
	
	# differentiate wrt u dfdu
	dfdu <- matrix(NA,x,m)
	f0 <- funx(X,u0=u0, PM,n,m)
	    for (i in 1:m)
	    {
	    xmi <- u0+Vu[,i]*dx
	    fi <- funx(X,u0=xmi, PM,n,m)
	    dfdu[,i]  <- spm_dfdx(fi,f0,dx)
	    }
      
    #  Create M0
    M0 <- rbind(rep(0,(x)+1),cbind(f0-dfdx%*%X,dfdx))
    #  Create M1
    M1<-rep(0,(length(M0)*m))
	dim(M1)<-c(dim(M0),m)
	for (i in 1:m){
	M1[,,i] <- rbind(rep(0,(x)+1),cbind((dfdu[,i]-dfdxu[,,i]%*%X),dfdxu[,,i]))
	}
    # Create L1
	dldx <- matrix(0,n,x)
	f0 <- funl2(X, PM,n,TE)
	    for (i in 1:x)
	    {
	    xmi <- X+Vx[,i]*dx
	    fi <- funl2(xmi, PM,n,TE)
	    dldx[,i]  <- spm_dfdx(fi,f0,dx)
	    }
      
    L1 <- cbind(l0-dldx%*%X,dldx)
    bireduce<-list()
    bireduce$M0 <- M0
    bireduce$M1 <- M1
    bireduce$L1 <- L1
    bireduce
    }

