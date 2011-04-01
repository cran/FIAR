dcmEstimate <-
function(DCM=DCM,ts){
    DCM<-stimfun(DCM)
    DCM$priors<-spm_dcm_priors(DCM) 
    pC <- DCM$priors$pC
    pE <- DCM$priors$pE 
    if(length(DCM$X0)==0){
#     X0 <- HPF(DCM)
#   
#     x0 <- matrix(rep(1,DCM$v))
#     if (length(X0)>1){ 
#     X0 <- as.matrix(cbind(x0,X0))
#     } else { X0 <- x0}}else{
    X0 <- matrix(rep(1,DCM$v))}else{
X0<-as.matrix(DCM$X0)}
    ## spm_nlsi
    nr <- length(ts)/DCM$v
    
    nh <- DCM$n
    nt <- length(ts)
    nq <- nr*DCM$v/nt
    h <- rep(0,nh)
    nb <- dim(X0)[1]
    nx <- nr*DCM$v/nb
    dfdu <- diag(nx)%x%X0
    Q <- diag(nx)%x%rep(1,DCM$v)
    hE <- rep(0,nh)
    ihC <- diag(nh)
    ##   ## svd priors$pC in R or ...##
    V <- spm_svd(pC,exp(-16))
    V2 <- Matrix(V)
    ## ... Copy decomposed matrix from Matlab
    
    ## spm_nlsi ##
    nu <- dim(dfdu)[2]
    np <- dim(V)[2]
    ip <- 1:np
    iu <- (1:nu)+np
    ## 2nd order moments
    pC    <- t(V2)%*%pC%*%V2
    uC    <- diag(nu)/1e-8
    ipC <- solve(bdiag(pC,uC))
    # initialize conditional density
    Eu    <- solve(t(dfdu)%*%dfdu)%*%(t(dfdu)%*%c(ts))
    p     <- rBind(t(V)%*%(pE-pE), Eu)
    Ep    <- pE + V%*%p[ip]
    
    Cp    <- pC
    
    ## EM
    CF   <- -Inf
    t    <- 256
     dFdh <- rep(0,nh)
     dFdhh <- matrix(0,nh,nh)
     Pr <- vector('list',length=nh)
     PS <- vector('list',length=nh)
    dx <- exp(-8)
    dfdp<-matrix(0,DCM$v*nr,np)
 options(warn=-1)
 for(ka in 1:64){# 1 ipv 128 om te testen


    ####  spm_int ##########
    
    f0 <- spm_int(PM=Ep,DCM=DCM)
    
    
	
	for (i in 1:dim(dfdp)[2]){
	xmi <- Ep+V2[,i]*dx
	fi <- spm_int(PM=xmi,DCM)
	dfdp[,i]  <- spm_dfdx(fi,f0,dx)
	}
    
    
    # prediction error and full gradients
    e     <- Matrix(c(ts) - c(f0) - dfdu%*%p[iu])
      
#	J     <- cBind(-dfdp, -dfdu)
    J     <- Matrix(cBind(-dfdp, -dfdu))
    
    
    ## M-STEP
    
	for (im in 1:16){
	#precision and conditional covariance

        iS <- Diagonal(nt)*1e-8
	    for (i in 1:nh){
	      iS <- iS + Diagonal(x=Q[,i])*exp(h[i])
	      }
	  S     <- solve(iS)
		iS    <- diag(nq)%x%iS
		Cp    <- solve(t(J)%*%iS%*%J + ipC) 
	        
      
   for (i in 1:nh){
	    Pr[[i]] <- Diagonal(x=Q[,i])*exp(h[i])
	    PS[[i]] <- Pr[[i]]%*%S
	    Pr[[i]] <- diag(nq)%x%Pr[[i]]
	    }
	
	# derivatives: dLdh = dL/dh,...
	   
    for (i in 1:nh){
	    dFdh[i] <- as.numeric(sum(diag(PS[[i]]))*nq/2-as.numeric(t(e)%*%Pr[[i]])%*%e/2-sum(Cp*(t(J)%*%Pr[[i]])%*%J)/2)
	    
		for (j in i:nh){
		dFdhh[i,j] <- -sum(colSums(as.matrix(PS[[i]]*PS[[j]])))*nq/2
		dFdhh[j,i] <-  dFdhh[i,j]
		}
	    }

    

	# add hyperpriors
	d     <- h - hE;
	dFdh  <- dFdh  - ihC%*%d;
	dFdhh <- dFdhh - ihC;  
	#update ReML estimate
	Ch    <- solve(-dFdhh)
	dh    <- Ch%*%dFdh
	h     <- h  + dh
	    for(i in 1:nh){
	    h[i]     <- max(h[i],-16)
	    }
	dF    <- t(dFdh)%*%dh
	if( as.numeric(dF) < 10^-2){break}
	}
    
    ## E-STEP
    
	
    F <- as.numeric(-1*t(e)%*%iS%*%e/2- crossprod(p,ipC)%*%p/2 - crossprod(d,ihC)%*%d/2 - DCM$v*nr*log(8*atan(1))/2 - logdet(S)*nq/2 + logdet(ipC%*%Cp)/2 + logdet(ihC%*%Ch)/2)
	
	if(F>CF){
	C_p   <- p
	C_h   <- h
	CF   <- F
	      # E-Step: Conditional update of gradients and curvature
	dFdp  <- t(-J)%*%iS%*%e - ipC%*%p
	dFdpp <- t(-J)%*%iS%*%J - ipC
	# decrease regularization
	t     <- t*2
	str <- paste("EM-step(-):")
	}else{
	p     <- C_p
	h     <- C_h
	    # and increase regularization
	t     <- min(t/2,128)
	str <- paste("EM-step(+):")
	}
	
    ## E-STEP update
    dp    <- spm_dx(dFdpp,dFdp,t)
    
    p  <- p + dp
    Ep <- pE + V%*%p[ip]

    ## Convergence
    dF <- as.numeric(crossprod(dFdp,dp))
    cat(str, ka,'   F:',CF,'  dF:',dF,"\n")
    if (ka > 2 && dF < 1e-2){break}
    }

    
    DCM$Ep <- Ep
	Ep <- Ep[-1]
	j <- 1:(DCM$n*DCM$n)
    DCM$A <- matrix(Ep[j],DCM$n,DCM$n)
	if(length(DCM$names)==DCM$n) {colnames(DCM$A)<-DCM$names
	rownames(DCM$A)<-DCM$names}
	
	Ep <- Ep[-j]
    j <- 1:(DCM$n*DCM$n*DCM$m)
    B0 <- Ep[j]
    dim(B0) <- c(DCM$n,DCM$n,DCM$m)
    DCM$B <- B0
if(length(DCM$names)==DCM$n){
for (z in 1:dim(DCM$B)[3]){
colnames(DCM$B[,,z])<-DCM$names
	rownames(DCM$B[,,z])<-DCM$names}
}
    Ep <- Ep[-j]
    j <- 1:(DCM$n*DCM$m)
    DCM$C <- t(matrix(Ep[j],DCM$n,DCM$m))
    if(length(DCM$names)==DCM$n) {
    colnames(DCM$C)<-DCM$names
    rownames(DCM$C) <- names(DCM$ons)}
    Ep <- Ep[-j]
    DCM$H <- matrix(Ep,nrow=DCM$n)   
DCM$Cp <- as.matrix(V%*%Cp[ip,ip]%*%t(V))
    PP <- Ncdf(DCM$Ep,diag(DCM$Cp))
 
	PP <- PP[-1]
	j <- 1:(DCM$n*DCM$n)
    DCM$pA <- matrix(PP[j],DCM$n,DCM$n)
	if(length(DCM$names)==DCM$n) {colnames(DCM$pA)<-DCM$names
	rownames(DCM$pA)<-DCM$names}
	
	PP <- PP[-j]
    j <- 1:(DCM$n*DCM$n*DCM$m)
    BB <- PP[j]
    dim(BB) <- c(DCM$n,DCM$n,DCM$m)
    DCM$pB <- BB
if(length(DCM$names)==DCM$n){
for (z in 1:dim(DCM$pB)[3]){
colnames(DCM$pB[,,z])<-DCM$names
	rownames(DCM$pB[,,z])<-DCM$names}
}

    PP <- PP[-j]
    j <- 1:(DCM$n*DCM$m)
    DCM$pC <- t(matrix(PP[j],DCM$n,DCM$m))
if(length(DCM$names)==DCM$n) {
    colnames(DCM$pC)<-DCM$names
    rownames(DCM$pC) <- names(DCM$ons)}
    DCM$F <- F
    DCM$Ce <- as.matrix(S)
    DCM
   
    }

