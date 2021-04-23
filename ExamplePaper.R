


source("lognormalSumCMC.RCMC.R")

# 1. Simple examples how to us CMC() and CMC.RIS to estimate right and left tail probabilities

########################################
# Example how to use CMC() CMC.RIS() Left Tail probabilities
CMC(n=1.e5,gam=5,d=10,sigma=1,rho=0.2)
CMC.RIS(n=1.e4,gam=5,d=10,rho= 0.2,sigma=1)
CMC.RIS(n=1.e4,gam=5,d=10,rho= 0.2,sigma=1,nin=10)# uses a larger number of RIS repetitions




##########################################
# 2. Comparison example of the paper for the left tail:
#    for 18 parameter settings d= c(4,10,30), s=c(0.25,1), rho=c(-0.5/d-1,0,0.5) 
#    with 3 gamma (ie. threshold) for each parameter setting such that the
#    left tail probabilities are approximately 1.e-3, 1.e-10 and 1.e-30 

load(file="Ex1.parm.Rdata")) 
#loads the 54 x 4 matrix Ex1.parm holding 54 rows with d, s, rho, gamma values


ManyExperGam.m <- function(simFun=CMC,pm=Ex1.parm,n=1.e4,...){
# makes for simulation method simFun the experiments with the parameters in pm
# for sample size n.
# returns a data.frame holding the input parameters and the simulated values:
#  d    s  rho gamma          est           SE      relErr seconds (ie used CPU time)
#  4 0.25 -0.5 3.109 1.007645e-03 1.792579e-06 0.001778979    0.06  
resm <- NULL
for( i in 1:length(pm$d)){
#    	print(paste("d=",d=pm$d[i],"s=",pm$s[i],"rho=",pm$rho[i]))
		t.0<-proc.time()[1]
		res <- simFun(n=n,d=pm$d[i],sigma=pm$s[i],rho=ifelse(pm$rho[i]<0,pm$rho[i]/(pm$d[i]-1),pm$rho[i]),gamma = pm$gamma[i],...) 
		resm<-rbind(resm,c(res,seconds= unname(proc.time()[1]-t.0))) 
}
data.frame(cbind(pm,resm))
}
ManyExperGam.m(simFun=CMC.RIS,pm=Ex1.parm[1:9,],n=1.e3) # 9 experiments with n only 1000

# comparison experiments:
system.time(resCMC <-ManyExperGam.m(simFun=CMC,pm=Ex1.parm,n=1.e5)) # last about 10 seconds ( for 54 CDF estimates)
#   user  system elapsed 
#   8.75    0.38    9.14 
system.time(resCMC.RIS <-ManyExperGam.m(simFun=CMC.RIS,pm=Ex1.parm,n=1.e4,nin=4))# last about 10 seconds ( for 54 CDF estimates)
#   user  system elapsed 
#   6.49    0.18    6.65 

# comparison of relative errors:

c(min(resCMC$relErr),median(resCMC$relErr),max(resCMC$relErr))
# [1] 0.0002320282 0.0032530798 0.2248736975
c(min(resCMC.RIS$relErr),median(resCMC.RIS$relErr),max(resCMC.RIS$relErr))
# [1] 0.0004258808 0.0012257987 0.0082423014

# Comparison of WNRelErr (Work Normalized Relative Error
WNRelErr<-resCMC$relErr*sqrt(resCMC$seconds) 
c(min(WNRelErr),median(WNRelErr),max(WNRelErr)) 
#[1] 0.0001030585 0.0010021826 0.1917518312   # WNRelErr for CMC
WNRelErr<-resCMC.RIS$relErr*sqrt(resCMC.RIS$seconds) 
c(min(WNRelErr),median(WNRelErr),max(WNRelErr))
#[1] 0.0001809518 0.0003510701 0.0038021145 # WNRelErr for CMC.RIS nin=4


################################
# Experiments to compare the relative errors of CDF and pdf estimates.

ManyExperCDFpdf.m <- function(simFun=CMC,pm=Ex1.parm,n=1.e4,...){
# makes for simulation method simFun the experiments with the parameters in pm
# for sample size n.
# returns a data.frame holding the input parameters and the simulated values of CDF and pdf:
CDFm <- pdfm <- NULL
for( i in 1:length(pm$d)){
#    	print(paste("d=",d=pm$d[i],"s=",pm$s[i],"rho=",pm$rho[i]))
		t.0<-proc.time()[1]
		res <- simFun(n=n,d=pm$d[i],sigma=pm$s[i],rho=ifelse(pm$rho[i]<0,pm$rho[i]/(pm$d[i]-1),pm$rho[i]),gamma = pm$gamma[i],...) 
		CDFm<-rbind(CDFm,c(res[1:3])) 
		pdfm<-rbind(pdfm,c(res[4:6],seconds= unname(proc.time()[1]-t.0))) 
}
data.frame(cbind(CDFm,pdfm))
}
# Experiment to demonstrate the behaviour of the PDF-estimates of CMC.RIS:
res<-ManyExperCDFpdf.m(simFun=CMC.RIS,pm=Ex1.parm,n=1.e5,pdfyn=T)
# to checks if the relative errors of pdf and CDF estimates are similar:
max(res$relErr/res$reEr.pdf)
#[1] 1.056486
min(res$relErr/res$reEr.pdf)
#[1] 0.9659838


# Experiment to demonstrate the behaviour of the PDF-estimates of CMC:
res<-ManyExperCDFpdf.m(simFun=CMC,pm=Ex1.parm,n=1.e6,pdfyn=T)
# to check if the relative errors of pdf and CDF estimates are similar:
max(res$relErr/res$reEr.pdf)
#[1] 1.182084
min(res$relErr/res$reEr.pdf)
#[1] 1.005113

 



