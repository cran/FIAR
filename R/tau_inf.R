tau_inf <-
function(data,p=.05,bin=0){

l <- ncol(data)
r <- nrow(data)
out<-list()
out$prb <- matrix(0,l,l)
out$orig <- matrix(0,l,l)

for(ii in 1:(l-1)) {
for (jj in (ii+1):l){

out$orig[ii,jj]=tau(data[,ii],data[,jj],bin=bin)
out$orig[jj,ii]=-out$orig[ii,jj]
}}

for(ii in 1:(l-1)) {
for (jj in (ii+1):l){
if (quantile(t,1-p/(ncol(data)^2-ncol(data)))>out$orig[ii,jj])
{out$prb[jj,ii]=0} else {out$prb[jj,ii]=1}

if (quantile(t,1-p/(ncol(data)^2-ncol(data)))>out$orig[jj,ii])
{out$prb[ii,jj]=0} else {out$prb[ii,jj]=1}

}}
out$orig<-t(out$orig)
out
}
