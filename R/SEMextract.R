#' @export
SEMextract<-function(ts,ons,dur,TR){
semts <- numeric()
for (l in 1:length(ons)){
time <- (ons[l]:(ons[l]+dur))
semts <- c(semts,time)
}
tsf <- ts[sort(semts+6/TR)[which(sort(semts+6/TR)<=nrow(ts))],]
tsf
}

