funl2 <-function(x,P,n,TE){
    x1 <- matrix(x,nrow=n)
    H <- P[(length(P)-n*6+1):length(P)]
    H0 <- matrix(H,nrow=n)
    V0        <- 100*0.04
    # slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
    # saturation Y:  R_iv = r0*[(1-Y)-(1-Y0)]
    r0        <- 25 # [Hz]
    # frequency offset at the outer surface of magnetized vessels
    nu0      <- 40.3 #% [Hz]
    # Get estimates of hemodynamic parameters
    E0       <- H0[,5]
    # estimated region-specific ratios of intra- to extravascular components of
    # the gradient echo signal (prior mean = 1, log-normally distributed scaling factor) 
    epsilon   <- exp(H0[,6])
    #% coefficients in BOLD signal model
    k1       <- 4.3*nu0*E0*TE
    k2      <- epsilon*r0*E0*TE
    k3      <- 1 - epsilon
    #% exponentiation of hemodynamic state variables
    x1[,2:5] <- exp(x1[,2:5]) 
    # output equation of BOLD signal model
    v <- x1[,4];
    q        <- x1[,5];
    y        <- V0*(k1*(1-q) + k2*(1-(q/v)) + k3*(1-v));
    y
    }

