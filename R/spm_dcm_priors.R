
spm_dcm_priors <-
function(DCM){
    a0 <- ifelse (DCM$a!=0,1,0)
    b0 <- ifelse (DCM$b!=0,1,0)
    c0 <- ifelse (DCM$c!=0,1,0)
    p <- 1e-3
    b <- 1
    s <- 3.0902 #fixed, but never used?
    q <- qchisq((1-p),DCM$n*(DCM$n-1))
    q     <- DCM$n/((DCM$n - 1)*q)
    #Intrinsic connections
    a0 <- matrix(a0,ncol=DCM$n,byrow=TRUE)
    diag(a0) <- 0
    dim(b0) <- c(DCM$n,DCM$n,DCM$m)
for (i in 1:DCM$m){
    b0[,,i] <- t(b0[,,i])*q
              }
    pC <- diag(c(1/16, a0*q, b0, c0)) # Covariance neuronal state
    #Expectations
    Ae <- diag(-1,DCM$n)
    Be <- ifelse (DCM$b!=0,0,0)
    Ce <- ifelse (DCM$c!=0,0,0)
    pE <- cbind(c(log(b),Ae,Be,Ce)) # Expectation neuron state

    qE <- DCM$h # Expectation hemod. state
    qC<- matrix(c(1.4837866e-2,-1.2823236e-4,-4.0541540e-4,3.4889257e-4,-3.4493093e-4,0,
    -1.2823236e-4,1.8370175e-3,-7.9879877e-4,4.2061309e-4,-1.1577819e-5,0,
    -4.0541540e-4,-7.9879877e-4,5.1387416e-2,2.0607856e-3,1.7781387e-3,0,
    3.4889257e-4,4.2061309e-4,2.0607856e-3,2.0564973e-4,6.6429753e-5,0,
    -3.4493093e-4,-1.1577819e-5,1.7781387e-3,6.6429753e-5,6.9002931e-5,0,
    0,0,0,0,0,3.1250000e-2),ncol=6)
    qC <- qC%x%diag(1,DCM$n)
    qE <- qE%x%cbind(rep(1,DCM$n))
    pE <- rbind(pE,qE)

    pC <- as.matrix(bdiag(pC,qC))
    p<-list()
    p$pC<-pC
    p$pE<-pE
    p$qC<-qC
    p$qE<-qE
    p
    }
