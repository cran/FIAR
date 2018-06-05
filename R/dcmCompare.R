#' @export
dcmCompare <-
function(DCM1, DCM2){

    bf<-list()
    cat('\n','AIC overall Bayes Factor','\n\n')
    nats <-  -1*(DCM1$AIC-DCM2$AIC)
    bits <- nats/log(2)
    bf_aic <- 2^(-bits)
    bf$bf_aic <- bf_aic
    cat('      BF: ',bf_aic,'\n')

    cat ('\n','BIC overall Bayes Factor','\n\n')
    nats <-  -1*(DCM1$BIC-DCM2$BIC)
    bits <- nats/log(2)
    bf_bic <- 2^(-bits)
    bf$bf_bic <- bf_bic
    cat('      BF: ',bf_bic,'\n\n')
    #bf
    }

