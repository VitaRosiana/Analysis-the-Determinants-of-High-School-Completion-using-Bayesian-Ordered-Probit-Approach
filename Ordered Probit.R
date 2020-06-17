
#--------------Ordered (Probit binary)----------------------------#
rm(list=ls())
#setwd("H:/STAT7016/Final Project/Final Project")
#setwd("C:/SEMESTER 4/STAT7016/Final Project")

#jae.4486.5.data.2 <- read.table("C:/SEMESTER 4/STAT7016/Final Project/jae-4486-5-data-2.txt", quote="\"", comment.char="")
#jae.4486.5.data.2 <- read.table("H:/STAT7016/Final Project/Final Project/jae-4486-5-data-2.txt", quote="\"", comment.char="")
jae.4486.5.data.2 <- read.table("C:/SEMESTER 4/Final Project/Final Project/jae-4486-5-data-2.txt", quote="\"", comment.char="")
data.probit <- jae.4486.5.data.2 
names(data.probit) <- c("school_ID", "parental_income", "cog_test","fathers_edu", "unknown_fathedu", "mothers_edu",
                        "unknown_motheredu", "num_siblings", "unknown_num_siblings", "female", "minority", "employment_rate",
                        "age_on_january","elig_dropout", "postsecond_edu", "grade_completion", "unemployed_stat")
names(data.probit)
#y <- data.probit[,17]
#y.new <- ifelse(y > 0.5,1,0)
completion <- rep(1,length(data.probit[,16]))
completion <- ifelse(data.probit[,16]==0,4,completion)
completion <- ifelse(data.probit[,16]==1,1,completion)
completion <- ifelse(data.probit[,16]==2,2,completion)
completion <- ifelse(data.probit[,16]==3,3,completion)
completion
y <- completion
#X <- data.probit[,2:11]
#X <- data.probit[,2:5]
X <- data.probit[, c(2,3,4,6,8,10,11)] # 7 covariates
X <- as.matrix(X)

#do.ninth <-ifelse(data.probit[,16]==1,1,0)
#do.tenth<-ifelse(data.probit[,16]==2,1,0)
#do.eleventh<-ifelse(data.probit[,16]==3,1,0)
#complete <- ifelse(data.probit[,16]==4)
#X<-cbind(rep(1,length(y)),X,do.ninth, do.tenth,do.eleventh)
#X <- as.matrix(X)

install.packages("MCMCpack")
library(MCMCpack)
deg.mcmc<-MCMCoprobit(y ~ X,mcmc=25000)
plot(deg.mcmc)
summary(deg.mcmc)


#-------------Ordered Probit Codes from the Lecture Notes ----------------#
rm(list=ls())
setwd("H:/STAT7016/Final Project/Final Project")
#setwd("C:/SEMESTER 4/STAT7016/Final Project")
#jae.4486.5.data.2 <- read.table("C:/SEMESTER 4/STAT7016/Final Project/jae-4486-5-data-2.txt", quote="\"", comment.char="")
jae.4486.5.data.2 <- read.table("H:/STAT7016/Final Project/Final Project/jae-4486-5-data-2.txt", quote="\"", comment.char="")
data.probit <- jae.4486.5.data.2 
names(data.probit) <- c("school_ID", "parental_income", "cog_test","fathers_edu", "unknown_fathedu", "mothers_edu",
                        "unknown_motheredu", "num_siblings", "unknown_num_siblings", "female", "minority", "employment_rate",
                        "age_on_january","elig_dropout", "postsecond_edu", "grade_completion", "unemployed_stat")
names(data.probit)
#y <- data.probit[,17]
#y.new <- ifelse(y > 0.5,1,0)
completion <- rep(1,length(data.probit[,16]))
completion <- ifelse(data.probit[,16]==0,4,completion)
completion <- ifelse(data.probit[,16]==1,1,completion)
completion <- ifelse(data.probit[,16]==2,2,completion)
completion <- ifelse(data.probit[,16]==3,3,completion)
completion
y <- completion
#X <- data.probit[,2:11]
#X <- data.probit[,2:5]
X <- data.probit[, c(2,3,4,6,8,10,11)] # 7 covariates
X <- as.matrix(X)

#X <-cbind(ychild,ypdeg,ychild*ypdeg)
X <- X
#y<-ydegr
y <- y
keep<- (1:length(y))[ !is.na( apply( cbind(X,y),1,mean) ) ]
X<-X[keep,] ; y<-y[keep]
#ranks<-match(y,sort(unique(y))) ; 
ranks<-y
uranks<-sort(unique(ranks))
K<-length(uranks)
n<-dim(X)[1] ; p<-dim(X)[2]
iXX<-solve(t(X)%*%X)  ; V<-iXX*(n/(n+1)) ; cholV<-chol(V)

#assumed prior on g's
mu<-rep(0,K-1) ; sigma<-rep(100,K-1)

###starting values
set.seed(1)
beta<-rep(0,p) 
z<-qnorm(rank(y,ties.method="random")/(n+1)) # divide by n+1 so that ranks are between 0 and 1.
g<-rep(NA,length(uranks)-1) #K-1 thresholds if there are K unique ranks
K<-length(uranks)



#objects to store results
BETA<-matrix(NA,20000,p) ; Z<-matrix(NA,20000,n) ; ac<-0

S<-40000
for(s in 1:S) 
{
  
  #update g 
  for(k in 1:(K-1)) 
  {
    a<-max(z[y==k]) #constraints on g_k
    b<-min(z[y==k+1])
    u<-runif(1, pnorm( (a-mu[k])/sigma[k] ),  #inverse CDF methods to draw a random g_k
             pnorm( (b-mu[k])/sigma[k] ) )
    g[k]<- mu[k] + sigma[k]*qnorm(u)
  }
  
  #update beta
  E<- V%*%( t(X)%*%z )
  beta<- cholV%*%rnorm(p) + E
  
  #update z
  ez<-X%*%beta
  a<-c(-Inf,g)[ match( y-1, 0:K) ] #match returns a vector of the positions of (first) matches of its first argument in its second.
  #1 is matched with zero; 2 is matched with 1.....(defines lower limit)
  b<-c(g,Inf)[y]  #assign upper limit to each z_i (i=1,..n)
  u<-runif(n, pnorm(a-ez),pnorm(b-ez) )
  z<- ez + qnorm(u)
  
  
  #help mixing
  c<-rnorm(1,0,n^(-1/3))  
  zp<-z+c ; gp<-g+c #add extra random component to zp and gp
  lhr<-  sum(dnorm(zp,ez,1,log=T) - dnorm(z,ez,1,log=T) ) + 
    sum(dnorm(gp,mu,sigma,log=T) - dnorm(g,mu,sigma,log=T) )
  if(log(runif(1))<lhr) { z<-zp ; g<-gp ; ac<-ac+1 }
  
  if(s%%(S/20000)==0) 
  { 
    cat(s/S,ac/s,"\n")
    BETA[s/(S/20000),]<-  beta
    Z[s/(S/20000),]<- z
  }
  
} 

par(mfrow=c(2,3))
blabs<-c(expression(beta[1]),expression(beta[2]),expression(beta[3]),expression(beta[4]), expression(beta[5]),
         expression(beta[6]), expression(beta[7]))
#thin<-c(1,(1:1000)*(S/1000))
#thin<-c(1,(1:250)*(10000/250))
#thin<-c(1,(1:200)*(10000/200))
thin<-c(1,(1:1000)*(20000/1000))

#Diagnostic check for posterior beta_j
for (j in 1 : p ){
  plot(thin,BETA[thin,j],type="l",xlab="iteration",ylab=blabs[j])
  abline(h=mean(BETA[,j]) )  
}

par(mfrow=c(2,3))
for (i in 1:p){
  acf(BETA[thin,i],xlab="iteration",ylab=blabs[j] )
}

beta.post <- apply(BETA, 2, mean)
beta.post

ZPM<-apply(Z,2,mean)
ZPM

# Posterior predicitive check
y.post <- rep(1, length(y))
y.post <- ifelse(ZPM < g[1] & ZPM > Inf, 1, y.post)
y.post <- ifelse(ZPM < g[2] & ZPM > g[1], 2, y.post)
y.post <- ifelse(ZPM < g[3] & ZPM > g[2], 3, y.post)
y.post <- ifelse(ZPM < Inf & ZPM > g[3], 4, y.post)

#Misclassification measurement
1-(sum(y.post==y)/length(y))

# Confidence Interval
apply(BETA,2,function(x) quantile(x,prob=c(.025,.5,.975)))

#mcmcpack
deg.mcmc<-MCMCoprobit(y ~ X,mcmc=25000)
plot(deg.mcmc)
summary(deg.mcmc)

