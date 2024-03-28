
library('lsei')
aa=c('A','C','D','E','F','G','H','I','K','L',
     'M','N','P','Q','R','S','T','V','W','Y')
# neh=c(4,1.62000000476837,1.88999998569489,3.95000004768372,0.319999992847443,
#       2.05999994277954,2.88000011444092,1,0.540000021457672,0.600000023841858,
#       1.12000000476837,1.88999998569489,2.58999991416931,3.95000004768372,
#       3.4300000667572,2.60999989509583,0.200000002980232,0.560000002384186,
#       0.0799999982118607,0.419999986886978)

neh=c(2, 3, 3, 2, 1, 4, 1, 2, 2, 3, 0, 1, 0, 2, 1, 1, 4, 0, 3, 1)

path="C:\\Workplace\\Python\\AnalysisForThePaper\\NEH\\Final\\Simulation_ape_mpe_asymp\\systmatic_error_check\\"

A=read.csv(paste(path,"A.csv",sep = ''))
# b=read.csv(paste(path,"b.csv",sep = ''))
b=read.csv(paste(path,"b_true.csv",sep = ''))

A=as.matrix(A)
b=as.matrix(b)
temp_res=list()
df2=data.frame(aa,neh)
for (i in seq(1:1000)){
  
  b=read.csv(paste(path,"b_true.csv",sep = ''))
  b=as.matrix(b)
  
  # error= rnorm(dim(A)[1],mean = 0,sd = 3)
  error= rnorm(dim(A)[1],mean = 1,sd = 3)
  b=b+error
  
  # res=solve(t(A)%*%A ) %*% t(A) %*% b 
  res=nnls(A,b)$x
  df2[paste("test ",i)]<- res;#list(100*abs(neh-res)/res)  
  
}

# write.csv(df2, paste(path,"res_with_no_noise.csv",sep=""), row.names = FALSE)
write.csv(df2, paste(path,"res_with_shift.csv",sep=""), row.names = FALSE)
# write.csv(df2, paste(path,"res_with_no_shift.csv",sep=""), row.names = FALSE)

