# Test.CMC
Efficient R Code for the CDF of the sum of (iid and exchangeable) lognormal random variables.
The code includes two functions:
CMC(n=1.e5,rep.out=1,gamma=0.4,d=4,sigma=1,rho=0.1,lower.tail=T) # simple and fast

CMC.RIS(n=1.e4,rep.out=1,d=11,rho= 0,sigma=1,gamma=d*0.8,lower.tail=T,nStep=5,nin=4,etnormtheta=2)
# more complicated but much more precise for very small tail probabilities (below 1.e-3)

The Algorithms are decribed in all detail in the working paper:
Dingeç, K. D. and Hörmann, W
Efficient Algorithms for Tail Probabilities of Exchangeable Lognormal Sums
