
aa=c('A','C','D','E','F','G','H','I','K','L',
     'M','N','P','Q','R','S','T','V','W','Y')
neh=c(4.3,2.3,2.3,5,3.3,2.5,3.75,1.25,4.5,1.25,
      1.25,2.5,3.75,5,2.3,3.75,0,1.25,0,1.25)

A=read.csv("C:\\Workplace\\Python\\AnalysisForThePaper\\NEH\\Final\\Simulation_ape_mpe_asymp\\A.csv")
# b=read.csv("C:\\Workplace\\Python\\AnalysisForThePaper\\NEH\\Final\\Simulation_ape_mpe_asymp\\b.csv")
b=read.csv("C:\\Workplace\\Python\\AnalysisForThePaper\\NEH\\Final\\Simulation_ape_mpe_asymp\\b_true.csv")

A=as.matrix(A)
b=as.matrix(b)

error= rnorm(dim(A)[1],mean = -1,sd = 3)
# error= rnorm(dim(A)[1],mean = 0.22341388086111777,sd = 1.4179280540373156)
b=b+error

res=solve(t(A)%*%A ) %*% t(A) %*% b

# plot(b-error,b)
# lines(b,b,col='red')

df=data.frame(aa,neh,res,100*abs(neh-res)/neh)

plot(df$neh,pch = 16)
axis(side = 1,at=seq(1:20), labels = aa)
points(df$X0,col='red')
axis(side = 1,at=seq(1:20), labels = aa)

df
