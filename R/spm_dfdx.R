spm_dfdx <-
function(f,f0,dx){
    dfdx  <- (as.vector(f) - as.vector(f0))/dx
    dfdx
    }

