
NEH_Compute <- function(A0, A1, A0_0, A1_0, pX)
{
  pH = 1.5574*10^(-4)
  NEH = (1 - pH - pX)*(1 - pH) / pX;
  NEH = NEH * (A1 / A0 - A1_0 / A0_0);
  return(NEH)
}

I0_t <- function(I0_0, pX, NEH)
{
  pH = 1.5574*10^(-4)
  I0t = I0_0 * (1 - pX / (1 - pH))^NEH;
  return(I0t)
  
}

I1_t <- function(I0_0, I1_0, NEH, pX)
{
  pH = 1.5574*10^(-4)
  I1t = (1 - pX / (1 - pH))^(NEH - 1) * pX * NEH * I0_0 / (1 - pH)^2;
  I1t = I1t + I1_0 * (1 - pX / (1 - pH))^NEH ;
  return(I1t);
}

I2_t <- function(I0_0, I1_0, I2_0, pX, NEH)
{
  
  pH = 1.5574*10^(-4)
  
  I0t = I0_t(I0_0, pX, NEH);
  
  I1t = I1_t(I0_0, I1_0, NEH, pX);
  
  I2t = I0t * (I2_0 / I0_0 - I1_0 * pH * NEH / (I0_0 * (1 - pH)));
  
  I2t = I2t + I0t * (pH / (1 - pH))^2 * NEH * (NEH + 1) / 2;
  
  I2t = I2t - I0t * ((pX + pH)/(1 - pH - pX))^2 * NEH * (NEH + 1) / 2;
  
  I2t = I2t + NEH * (pX + pH) * I1t / (1 - pH - pX);
  
  return(I2t);
}

Compute_neh_pxt <- function(I0_0,I1_0,I2_0,I0,I1,I2){ 
  
  
  pH = 1.5574*10^(-4)  
  
  delta_A2_t = I2 / I0 - I2_0 / I0_0;
  
  delta_A1_t = I1 / I0 - I1_0 / I0_0;
  
  A = delta_A2_t - (I1_0 / I0_0) * pH * delta_A1_t - (I1 / I0) * (1 - pH) * delta_A1_t;
  
  A = A + (-delta_A1_t^2 * (2*pH - 1) + delta_A1_t * ( 2 * pH - 1) / (1 - pH) ) / 2;
  
  B = -delta_A2_t * (1 - pH) + (I1_0 / I0_0) * pH * delta_A1_t * 2 * ( 1 - pH);
  
  B = B + (1 - 2 * pH) * I1 / I0 * ( 1 - pH) * delta_A1_t;
  
  B = B + ( delta_A1_t^2 * (1 - pH) * ( 2* pH - 1) + delta_A1_t^2 * 2 * pH * (1 - pH) - delta_A1_t * 2 * pH ) / 2;
  
  new_pxt=-B/A
  new_neh=NEH_Compute(I0,I1, I0_0,I1_0,pX = new_pxt );
  # print(c(A,B,new_neh))
  
  return (c(new_pxt, new_neh,A,B));
  
}



data=read.csv("C:\\Workplace\\C#\\Test\\Test_test\\ConsoleApp2\\bin\\Debug\\0_04.csv")

noise_levels=c(0,0.0000005,0.000005,0.00005,0.0005)
for(i in seq(1,dim(data)[1])){
  # print(data$NEH[i])
  # vals=Compute_neh_pxt(data$M0[i],data$M1[i],data$M2[i], data$I0[i],data$I1[i],data$I2[i]);
  
  # print('with noise')
  # noise=rnorm(6,0,noise_level)
  I0_N=data$I0[i];
  I1_N=data$I1[i];
  I2_N=data$I2[i] - data$I2[i]*4/100;
  I3_N=data$I3[i];
  I4_N=data$I4[i];
  I5_N=data$I5[i];
  
  sum_is=sum(I0_N+I1_N+I2_N+I3_N+ I4_N+I5_N)
  I0_N=I0_N/sum_is
  I1_N=I1_N/sum_is
  I2_N=I2_N/sum_is
  I3_N=I3_N/sum_is
  I4_N=I4_N/sum_is
  I5_N=I5_N/sum_is
  
  vals2=Compute_neh_pxt(data$M0[i],data$M1[i],data$M2[i], I0_N,I1_N,I2_N);
  if (vals2[3]>0  || vals2[4]<0 || vals[1]<0 || vals2[1]> 0.1) next;
  print(c(data$NEH[i],vals2[2],vals2[1] ))

}