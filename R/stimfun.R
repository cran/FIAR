stimfun <-
function(DCM){
    TR<-DCM$TR
    T<-DCM$T
    v <- DCM$v
    sf<-list()
	for (i in 1:length(DCM$ons)){
    
    if(length(DCM$dur[[i]])==1)
         {DCM$dur[[i]]<-rep(DCM$dur[[i]],length(DCM$ons[[i]]))}
    dur<-DCM$dur[[i]]
    ons<-DCM$ons[[i]]
    
      ifelse(dur==0,(u_sf<-T/TR),(u_sf<-1))
	ton <- round(ons*T)+32
	tof <- round(dur*T) + ton + 1
	sf[[i]] <- rep(0,DCM$v*T+128)
	ton <- pmax(ton,1)
	tof <- pmax(tof,1)
	    for(j in 1:length(ton))
	    {
	      if (length(sf[[i]])>ton[j]) sf[[i]][ton[j]]<-sf[[i]][ton[j]]+u_sf
	      if (length(sf[[i]])>tof[j]) sf[[i]][tof[j]]<-sf[[i]][tof[j]]-u_sf
	    }
    sf[[i]] <- cumsum(sf[[i]])
    sf[[i]] <- sf[[i]][1:(DCM$v*T + 32)]
    }
    stimfunc <-matrix(NA,length(sf[[1]]),length(DCM$ons))
      	for(i in 1:length(DCM$ons)){
			
	stimfunc[,i] <- sf[[i]]
	
                                 }
	#prepare sf for DCM (drop first 32 timepoints)
	
	sf <- as.matrix(stimfunc[33:nrow(stimfunc),])			 
	s <- ceiling((1:v)*nrow(sf)/v) # output times
    if (ncol(sf)>1){
    
    t <- NA
	for (i in 1:ncol(sf)){
		temp <- c(1,which(diff(sf[,i])!=0)+1)
	t <- c(t,temp) # input times
                           }
	t <- sort(t[-1])[-which(diff(sort(t[-1]))==0)]
	              }
		      
    if (ncol(sf)==1){
    
              t=c(1,which(diff(sf)!=0)+1) # input times
                      }
	T0 <- sort(c(s,t))
    s <- sort(c(s,t),index.return=TRUE)$ix
    dt0 <- c(TR/T*diff(T0),0)

DCM$sf<-sf
DCM$s<-s
DCM$T0<-T0
DCM$dt0<-dt0
DCM
}

