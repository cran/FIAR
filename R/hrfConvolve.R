hrfConvolve <-
function (x=NULL, scans = NA, onsets = c(), durations = c(), rt = NA, 
    SNR=0, mean = FALSE, a1 = 6, a2 = 12, b1 = 0.9, b2 = 0.9, 
    cc = 0.35) 
{
   hrf <- function(x, a1, a2, b1, b2, c) {
        d1 <- a1 * b1
        d2 <- a2 * b2
        c1 <- (x/d1)^a1
        c2 <- c * (x/d2)^a2
        res <- c1 * exp(-(x - d1)/b1) - c2 * exp(-(x - d2)/b2)
        res
    }
    if (is.null(x)){ 
    numberofonsets <- length(onsets)
    if (length(durations) == 1) {
        durations <- rep(durations, numberofonsets)
    }
    
 	stimulus <- rep(0, scans)
    for (i in 1:numberofonsets) {
        for (j in onsets[i]:(onsets[i] + durations[i] - 1)) {
            stimulus[j] <- 1
        }
    }
    hrfnn <- convolve(stimulus, hrf(scans:1, a1, a2, b1/rt,b2/rt, cc))
    
	if(SNR>0){
	sdS <- sd(hrfnn)	
	noise=rnorm(scans,sd=sdS/SNR)
	#Zx <- x/sd(x)
	#sdN <- sdS/SNR
	hrfnn <- hrfnn + noise
    }
	else{hrfnn <- hrfnn}
	
	if (mean) {
        hrfnn - mean(hrfnn)
    }
    else {
        hrfnn
		}
    }
	else{ hrfnn <- convolve(x, hrf(length(x):1, a1, a2, b1,b2, cc))
		
	if(SNR>0){
	sdS <- sd(hrfnn)	
	noise <- rnorm(length(x),sd=sdS/SNR)
	
	hrfnn <- hrfnn + noise
    }
	else{hrfnn <- hrfnn}
	
		if (mean) {
        hrfnn - mean(hrfnn)
    }
    else {
        hrfnn
		}
    }
				
}

