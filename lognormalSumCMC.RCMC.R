

########################################################################################################################
########################################################################################################################
#
# start of code of our two methods, CMC and CMC.RIS, 
# for left tail of sum of exchangable normals: 

k0k1 <- function(d=50,rho=0,sigma=1){
 alpha <- ifelse(rho==0, -(d-2)/2, (1-sqrt(1-rho*(rho*(d-1)-d+2)))/rho )
 k0 <- sqrt(d/(1+rho*(d-1)))/sigma
 k1 <- sigma*(alpha-1)/sqrt(alpha^2+d-1)
 cbind(k0,k1)
}
#k0k1(d=50,rho=0,sigma=0.2)

########################################################
# Clean fast Implementation of CMC algorithm 
#######################################
calcCiv <- function(uv){
# calculates the coefficients civ necesary to calculate the root
# uses the fastest implementation with complexity O(n*d)
# returns a  "n x d" matrix holding in its rows the coefficient vector civ
#
# uv ... is the same as zm, a matrix of N(0,1) variates
#        (technically in R a "n x (d-1) - matrix)

n <- length(uv[,1])
d <- length(uv[1,])+1
hv1 <- c( sqrt(1/((1:(d-1))*(2:d))) , 0 )
hv2 <- c( 0, sqrt((1:(d-1))/(2:d)) )
civm <- matrix(0,nrow=n,ncol=d);
civm[,d] <- -hv2[d]*uv[,d-1]
cusu <- numeric(n)
for(i in (d-1):2 ){
  cusu <- cusu+hv1[i]*uv[,i]
  civm[,i] <- cusu -hv2[i]*uv[,i-1] 
}
civm[,1] <- cusu + hv1[1]*uv[,1]
civm
}
#calcCiv(uv=rbind(c(0,1),cbind(1,2)))

CMC <- function(n=1.e5,rep.out=1,gamma=0.4,d=4,sigma=1,rho=0.1,pdfyn=F){
# simulates the CDF F(gamma) = Prob(S < gamma) for S sum of exchangeable lognormals
# for pdfyn=F a vector holding the CDF-estimate its SE and the relative error is returned
# for pdfyn=T also the pdf estimate its SE and relative error are appended
# using pdfyn=T also pdf(gamma) is estimated
# n ... sample size
# rep.out ... outer repetetitions used to reach very large sample sizes
#     rep.out =1 or rep.out >= 100 is recommended for reliable error estimates
# gamma ... threshold for the probabilities
# d, sigma, rho ... parameter of the exchangeable lognormal vector
# pdfyn = F ... only CDF values are calcuated, 
#       = T ... CDF and PDF values are calculated and returned    

kk <- k0k1(d=d,rho=rho,sigma=sigma) 
res.m <- matrix(NA,nrow=rep.out,ncol=ifelse(!pdfyn,3,6))
for(i in 1:rep.out){
  zm1 <- matrix(rnorm((d-1)*n),nrow=n)
  civm <- calcCiv(uv=zm1)
  root <- kk[1]*log(gamma/rowSums(exp(-kk[2]*civm)))
  CDF <- pnorm( root)
  if(!pdfyn){ 
     res <- c(est=est<-mean(CDF),SE=SE<-sd(CDF)/sqrt(n),relErr=SE/est)
	 res.m[i,] <- res
  }else{
     pdf <- dnorm( root)*kk[1]/gamma
     res <- c(est=est<-mean(CDF),SE=SE<-sd(CDF)/sqrt(n),relErr=SE/est,
						est.pdf=est.pdf<-mean(pdf),SE.pdf=SE.pdf<-sd(pdf)/sqrt(n),reEr.pdf=SE.pdf/est.pdf)
     res.m[i,] <- res  
  }
}	
if(rep.out==1) return(res)
res.CDF <- c( est = est<-mean(res.m[,1]) , SE = SE <-sd(res.m[,1])/sqrt(rep.out), relErr = SE/est)
if(!pdfyn) return(res.CDF)  
return(  c(res.CDF, est.pdf = est.pdf<-mean(res.m[,4]) , SE.pdf = SE.pdf <- sd(res.m[,4])/sqrt(rep.out), reEr.pdf = SE.pdf/est.pdf)  )
}
#res<-CMC(n=1.e6,gam=5,d=10,sigma=1,rho=0,pdfyn=T)
# CMC(n=1.e3,rep.out=1.e2,gam=.2,d=10,sigma=1,rho=-0.5/9,pdfyn=T)

####################################################
rcivm <- function(d=5,n=2){
#generates random civ vectors and returns them as rows of a "n x d"-matrix 
  zm1 <- matrix(rnorm((d-1)*n),nrow=n)
  zm1 <- zm1/sqrt(rowSums(zm1^2))
  civm=calcCiv(uv=zm1)
  civm
}#rcivm(d=4,n=2)
####



###############################################################################
dchilog <- function(x,df){
# log-density of chi distribution
 (df-1)*log(x)-x^2/2 -(df/2-1)*log(2)-lgamma(df/2)
}
###############################################################################
ddnorm <- function(xv){
# derivative of density of standard normal
 (-xv*exp(-0.5*xv^2)/sqrt(2*pi)) 
}

###############################################################################
# exponential Tail normal distribution used for IS
etnormSetUp <- function(theta){
# setup that returns the constants necessary for qISetnorm()
# etnorm... exponential tail normal distribution 
# etnorm  standard normal in the center and exponential |z|>theta
# thus the log density of etnorm is a line for $z$ > theta
# the first derivative of the log density is smooth, 
# the second derivative has a discontinuity at +/- theta
 pnth <- pnorm(-theta)
 In <- sqrt(2*pi)*(1-2*pnth)
 Ie <- 2*exp(-theta^2/2)/theta
 pe <- 0.5*Ie/(In+Ie)
 pn <- In/(In+Ie)
 c<- 1/(2/theta*exp(-theta^2/2)+sqrt(2*pi)*(1-2*pnth))
 c(theta,pnth,pe,pn,c)
}

qISetnorm <- function(u,cov){
# quantiles of etnorm used for inversion method when etnorm is applied for Importance sampling
# returns q(u) as first column vector and its density value as second column vector
# u ... vector of probabilities
# cov ... vector of constants (produced by etnormSetUp(theta=2))
theta<-cov[1];pnth<-cov[2];pe<-cov[3];pn<-cov[4];c<-cov[5]
ivc <- which( u > pe & 1-u > pe )
xv <- fv <- numeric(length(u))
xc <- qnorm(pnth+(u[ivc]-pe)/(1-2*pe)*(1-2*pnth))
fc <- c*exp(-xc^2/2)
ut <- u[-ivc]
xt <- ifelse(ut < 0.5,-1,1)*(theta - log(pmin(ut,1-ut)/pe)/theta)
ft <- c*exp(0.5*theta^2-theta*abs(xt))
xv[ivc] <- xc
fv[ivc] <- fc
xv[-ivc] <- xt
fv[-ivc] <- ft
cbind(xv,fv)
}
#cov2 <- etnormSetUp(theta=2)
#res<-qISetnorm(u=(1:99)/100,cov =cov3.5)


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
CMC.RIS <- function(n=1.e2,rep.out=1,gamma=5.3,d=11,rho= 0,sigma=1,nStep=5,
                    nin=4,etnormtheta=2,pdfyn=F){
# estimates the CDF value F(gamma) = Prob(sumlognormal < gamma)
# for pdfyn=F a vector holding the CDF-estimate its SE and the relative error is returned
# for pdfyn=T also the pdf estimate its SE and relative error are appended
# Details:
# applies CMC.RIS using a normal-exponential tail IS density for Radius CMC
# for one random direction "nStep"+1 evaluations of the optimal IS density
# are used to find the mode; than for IS "nin" random variates from the 
# etnorm distribution (normal distribution with exponential tail for |x| >etnormtheta )
# are generated.
# nStep ... number of iterations for Newton method
# rep.out ... outer repetetitions used to reach very large sample sizes
#     rep.out =1 or rep.out >= 100 is recommended for reliable error estimates
# gamma ... threshold for the probabilities
# d, sigma, rho ... parameter of the exchangeable lognormal vector
# nin ... inner repetitions of RIS
# etnormtheta... default= 2 (so IS density is N(0,1) for abs(z)<2 and exponential in the tails
############################################################

CMC.RIS.sim <- function(kk=mykk,gam,civm=myciv,nStep=5,nin=4,stratyn=T,etnormtheta=2,logisticyn=F,pdfyn=F){
# uses Newton method with nStep steps to find the mode of optIS and then makes Radius IS
# returns a vector of n iid simulation results
##########
lSEcivm <- function(rv,civm,deriv=c(0,1,2)){
# rv ... vector of length n or of length 1
# civm ... nxd matrix with civ in every row
# calculates the log of the sum of exp(civ*r) 
# deriv = 0 ... return only vector of function values, 
# deriv = 1 ... returns matrix with function value as first column and derivative as second
# deriv = 2 ... returns three column-matrix holding function value, first derivative and second derivative
# rv=1:5;civm=rcivm(n=5,d=4);deriv=2

 if(length(rv)>1 & ( length(civm[,1])!= length(rv) ) ){
   print(paste("lSEcivm(): Error !!! length of rv unequal to row number of civm!! exiting!!"))
   return(NULL)
 } 
 
 exm <- exp(rv * rbind(civm))
 SE <- rowSums(exm)			# sum exp() for each r-value
 f <- log( SE)
 if(deriv==0) return(  cbind(f) )
 cSE <-  rowSums(exm <- exm*civm) 			# sum of ci*exp() 
 if(deriv==1) return( cbind(f,cSE/SE) )
 c2SE <- rowSums(exm*civm) 		# sum of ci^2*exp() 
 return( cbind(f,cSE/SE,c2SE/SE-(cSE/SE)^2) ) 
}

##################################################
optIScivmFast  <-  function(rv,kk,gam,civm,pdfyn=F){
# calculates the optimal IS density of radius IS
 d <- length(civm[1,])
 root <- kk[1]*(log(gam)-lSEcivm(-kk[2]*rv,-civm,deriv=0))
 if(!pdfyn)return(pnorm(root)*exp(dchilog(rv,d-1)))
 return(cbind(pnorm(root)*exp(dchilog(rv,d-1)),
                dnorm(root)*exp(dchilog(rv,d-1))*kk[1]/gam ) )
}

loptIScivm <-function(rv,kk,gam,civm){
# calculates the logarithm of the optimal IS density of radius IS
# together with the first 2 log-derivatives
# rv ... vector of radius values of length n or of length 1
# kk vector of k0k1
# civm ... n x d   matrix with civ in every rows
 d <- length(civm[1,])
 lsem <- lSEcivm(-kk[2]*rv,-civm,deriv=2)
 root <- kk[1]*(log(gam)-as.vector(lsem[,1]))
 droot <- kk[1]*kk[2]*as.vector(lsem[,2])
 ddroot <- -kk[1]*kk[2]^2*as.vector(lsem[,3])
 prob <- pnorm(root)
 dnroot <- dnorm(root)
 dlprob <- dnroot * droot / prob
 ddlprob <- (  ddnorm(root)*droot^2+dnroot*ddroot ) / prob - (dnroot * droot / prob)^2
 loIS <- log(prob) + dchilog(rv,d-1)  
 dloIS <- dlprob + (d-2)/rv -rv
 ddloIS <- ddlprob - (d-2)/rv^2 -1
 data.frame(loIS=loIS,dloIS=dloIS,ddloIS=ddloIS)
}



n<- length(civm[,1])
d <- length(civm[1,])
r <- sqrt(d-2)/2  # starting value of Newton method for left tail
  for(i in 1:nStep){
    ois <- loptIScivm(r,kk,gam,civm)
    step <- ois$dloIS/ois$ddloIS
    r <- r - ifelse(step<r, step, r/2) #x = x0 - f(x0)/f'(x0)
  }
 mode <- r 
 ois <- loptIScivm(mode,kk,gam,civm) 
 s_est <- sqrt(-1/ois$ddloIS)
 sumy <- 0
 covEtnorm <- etnormSetUp(theta=etnormtheta)
 for(i in 1:nin){
    zfm <- qISetnorm(u=(i-1+runif(n))/nin , cov =covEtnorm)
	xv <- mode + s_est*zfm[,1]# zfm[,1]... random variates of etnorm distribution 
    xv <- ifelse(xv <1.e-10,0,xv)
    sumy <- sumy+ optIScivmFast(xv,kk=kk,gam=gam,civm,pdfyn=pdfyn)/(zfm[,2]/s_est)
  }
 y <- sumy/nin
}

###################
# "main code of the function
 kk <- k0k1(d=d,rho=rho,sigma=sigma)
 res.m <- matrix(NA,nrow=rep.out,ncol=ifelse(!pdfyn,3,6))
 for(i in 1:rep.out){
   civm <- rcivm(d=d,n=n)
   y <- CMC.RIS.sim(kk,gam=gamma,civm,nStep=nStep,nin=nin,etnormtheta=etnormtheta,pdfyn=pdfyn)
#   res.m[i,] <- c(est<-mean(y),SE<-sd(y)/sqrt(n), SE/est)
  if(!pdfyn){ 
     res <- c(est=est<-mean(y),SE=SE<-sd(y)/sqrt(n),relErr=SE/est)
	 res.m[i,] <- res
  }else{
     CDF <- y[,1]
     pdf <- y[,2]
     res <- c(est=est<-mean(CDF),SE=SE<-sd(CDF)/sqrt(n),relErr=SE/est,
						est.pdf=est.pdf<-mean(pdf),SE.pdf=SE.pdf<-sd(pdf)/sqrt(n),reEr.pdf=SE.pdf/est.pdf)
     res.m[i,] <- res  
  }

 }
 if(rep.out==1) return(res)
 res.CDF <- c( est = est<-mean(res.m[,1]) , SE = SE <-sd(res.m[,1])/sqrt(rep.out), relErr = SE/est)
 if(!pdfyn) return(res.CDF)  
 return(  c(res.CDF, est.pdf = est.pdf<-mean(res.m[,4]) , SE.pdf = SE.pdf <- sd(res.m[,4])/sqrt(rep.out), reEr.pdf = SE.pdf/est.pdf)  )
}
##############################################

########################################
# Example how to use CMC() CMC.RIS()   only CDF values
# CMC(n=1.e5,gam=5,d=10,sigma=1,rho=0.2)
# CMC.RIS(n=1.e4,gam=5,d=10,rho= 0.2,sigma=1,nin=4)

########################################
# Example how to use CMC() CMC.RIS()   CDF and pdf values
# CMC(n=1.e5,gam=2,d=10,sigma=1,rho=0.2,pdfyn=T)
# CMC.RIS(n=1.e4,gam=2,d=10,rho= 0.2,sigma=1,nin=4,pdfyn=T)



