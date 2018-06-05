
HPF <-
function(DCM){
    if(DCM$HPF > 0) {
	    n <- floor( 2*(DCM$v*DCM$TR)/DCM$HPF + 1 )
	    K <- spm_dctmtx(DCM$v, n)[,-1]
        K
                 		}
    else {
	    K <- 1
	    K
 		}
    }

