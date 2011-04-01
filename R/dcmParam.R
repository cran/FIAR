dcmParam <-
function(a=NA, b=NA, c=NA, ons = list(), dur = list(), v=NA, n=NA, m=NA, TR=NA, h = c(0.65, 0.41, 0.98, 0.32, 0.34, 0), names=c(), TE = 0.04, T = 16, x = 5 * n, HPF=0, auto = FALSE) {
   if(auto == TRUE){
    DCM        <- list()
    DCM$names  <- list()
    DCM$inputs <- list()
    DCM$ons    <- list()
    DCM$dur    <- list()
    DCM$T<- 16 
    DCM$HPF <- 0
    DCM$h <- c(.65, .41, .98, .32, .34,0)    
            
       n <- readline("enter number of regions: ")
       DCM$n <- as.numeric(n)
       
       DCM$x <- 5*DCM$n
       
       
       for (i in 1:DCM$n){
       r <- paste('enter name of region', i,':',sep=' ')
       name <- readline(r)
       DCM$names[[i]] <- name
                      }
       
       TR <- readline("enter TR (in seconds): ")
       DCM$TR <- as.numeric(TR)
       
       v <- readline("enter number of scans: ")
       DCM$v <- as.numeric(v)
       
       TE <- readline("enter TE (in seconds): ")
       DCM$TE <- as.numeric(TE)
       
       m <- readline("enter number of inputs: ")
       DCM$m <- as.numeric(m)
       
      
       for (i in 1:DCM$m){
       
       r <- paste('enter name of input', i,': ',sep=' ')
       name <- readline(r)
       DCM$inputs[[i]] <- name
       
       r <- paste('enter onsets of', name, '(in scans, space separated): ',sep=' ')
       onsets <- readline(r)
       DCM$ons[[name]] <- as.numeric(scan(textConnection(onsets), what="character", sep=" ",quiet=TRUE) )

       
       r <- paste('enter duration of', name, '(in scans, space separateed): ',sep=' ')
       duration <- readline(r)
       DCM$dur[[name]] <- as.numeric(scan(textConnection(duration), what="character", sep=" ",quiet=TRUE) )
       
   
       }
       closeAllConnections()
       
       u <- unlist(DCM$inputs)
       n <- unlist(DCM$names)
       A <-diag(DCM$n)-diag(DCM$n)
       colnames(A) <- unlist(DCM$names)
       rownames(A) <- unlist(DCM$names)
       
       for (k in 1:length(n)){
       n2 <- n[-k] 
       for (j in 1:length(n2)){
       r   <- paste('anatomical connection from' , n[k], 'to',n2[j],'(scalar): ' ,sep=' ')
       con <- readline(r)
       A[which(n==n2[j]),which(n==n[k])] <- as.numeric(con)
       
                              }
			     }
       DCM$a <-c(t(A))		     
       
       B <-list()
       for (m in 1:DCM$m){
       B[[m]] <-diag(DCM$n)-diag(DCM$n)
       colnames(B[[m]]) <- unlist(DCM$names)
       rownames(B[[m]]) <- unlist(DCM$names)}
      		 
       for (i in 1:DCM$m){
       for (k in 1:length(n)){
       n2 <- n[-k] 
       for (j in 1:length(n2)){
       r   <- paste('functional influence of', toupper(u[i]), 'from' , n[k], 'to',n2[j],'(scalar): ' ,sep=' ')
       con <- readline(r)
       B[[i]][which(n==n2[j]),which(n==n[k])] <- as.numeric(con)
      
       
                              }
			     } 
		 
			 } 
	for (i in 1:DCM$m){
	B[[i]] <- t(B[[i]])
	                  }
	DCM$b <- unlist(B)
	
	C <-matrix(nrow=DCM$m,ncol=DCM$n)
	for (i in 1:DCM$m){
	for (k in 1:length(n)){ 
        r   <- paste('direct influence of', toupper(u[i]), 'to',n[k],'(scalar): ' ,sep=' ')
        con <- readline(r)
        C[i,k] <- as.numeric(con)
                               }
			  }
        DCM$c <- c(t(C))
	
		#DCM$X0 <- 1
		DCM$sf <- stimfun(DCM) 
	#	DCM$y<-spm_dcm_gen(DCM=DCM,SNR=0)
		DCM
		}
	else{
		DCM<-list()
		DCM$a<-a
		DCM$b<-b
		DCM$c<-c
		DCM$h<-h
		DCM$ons<-ons
		DCM$dur<-dur
		DCM$T<-T
		DCM$TR<-TR
		DCM$TE<-TE
		DCM$m<-m
		DCM$v<-v
		DCM$n<-n
		#DCM$HPF<-0
		DCM$x<-x
		DCM$names<-names
		
	#DCM$HPF <- 0
	DCM <- stimfun(DCM) 
	DCM
	}
	
	}

