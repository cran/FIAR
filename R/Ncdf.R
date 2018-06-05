
Ncdf<-function(u,v){
x=0
u=as.matrix(abs(u))
v=as.matrix(v)
ad=as.matrix(c(2,2,2))
rd=max(ad)
as1=c(1,1,rep(1,(rd-ad[1])))
as2=c(dim(u),rep(1,(rd-ad[2])))
as3=c(dim(v),rep(1,(rd-ad[3])))
as=rbind(as1,as2,as3)
rs=pmax(as1,as2,as3)
xa=rowMeans(as)>1
F=rep(0,max(rs))
md=v>0
Q=which(v>0)
if (xa[1]){Qx=Q} else {Qx=1}
if (xa[2]){Qu=Q} else {Qu=1}
if (xa[3]){Qv=Q} else {Qv=1}

X=(x[Qx]-u[Qu])/sqrt(2*v[Qv])
F[Q]=1-(0.5 + 0.5*(2*pnorm(X*sqrt(2))-1))
F
}
