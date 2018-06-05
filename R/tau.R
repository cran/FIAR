
tau <- function(d1,d2,bin=0){
XT=length(d1)
grotMIN=quantile(d1,.10,type=2)
grotMAX=quantile(d1,.90,type=2)
d1b=pmax(pmin((d1-grotMIN)/(grotMAX-grotMIN),1),0)

grotMIN=quantile(d2,.10,type=2)
grotMAX=quantile(d2,.90,type=2)
d2b=pmax(pmin((d2-grotMIN)/(grotMAX-grotMIN),1),0)
if (bin > 0){
d1b=(d1b>bin)
d2b=(d2b>bin)}

theta1=t(d1b)%*%d2b/XT
theta2=t(d1b)%*%(1-d2b)/XT
theta3=t(d2b)%*%(1-d1b)/XT
theta4=t(1-d1b)%*%(1-d2b)/XT

#EEE=(theta1+theta2)*(theta1+theta3)
#max_theta1=min(theta1+theta2,theta1+theta3)
#min_theta1=max(0,2*theta1+theta2+theta3-1)

  if(theta2>theta3){tau_12=1-(theta1+theta3)/(theta1+theta2)} else
  {tau_12=(theta1+theta2)/(theta1+theta3)-1}

tau_12

                      }

#write.table(d1,file='/home/bjorn/Matlab/d1.dat',col.names=F,row.names=F)
#write.table(d2,file='/home/bjorn/Matlab/d2.dat',col.names=F,row.names=F)
