#' @export
dcmEvidence <-
function(DCM,ts){
    #DCM$X0 <- 1
    DCM$priors<-spm_dcm_priors(DCM)
        if(length(DCM$X0)==0){
#     X0 <- HPF(DCM)
#
#     x0 <- matrix(rep(1,DCM$v))
#     if (length(X0)>1){
#     X0 <- as.matrix(cbind(x0,X0))
#     } else { X0 <- x0}}else{
    X0 <- matrix(rep(1,DCM$v))}else{
X0<-as.matrix(DCM$X0)}
    y <- spm_int(PM=DCM$Ep,DCM=DCM)
    R <- ts-y
    R <- R - X0%*%solve(t(X0)%*%X0)%*%(t(X0)%*%R)
    pCdiag <- diag(DCM$priors$pC)
    wsel <- which(pCdiag!=0)
    evidence_region_cost <- numeric(DCM$n)
	for (i in 1:DCM$n){
	lambda <-  as.matrix(DCM$Ce)[i*DCM$v,i*DCM$v]
	evidence_region_cost[i] <- -.5*DCM$v*log(lambda)
	evidence_region_cost[i] <-  evidence_region_cost[i] - ((.5*t(R[,i])*(1/lambda))%*%diag(DCM$v))%*%R[,i]
	}
    evidence_aic_penalty <- length(wsel)
    evidence_bic_penalty <- .5*length(wsel)*log(DCM$v)

    DCM$AIC <- sum(evidence_region_cost)- evidence_aic_penalty
    DCM$BIC <- sum(evidence_region_cost)- evidence_bic_penalty

    DCM
    }

