
data=read.csv("C:\\Workplace\\Python\\AnalysisForThePaper\\NEH\\second_round\\prepare_Data\\data.csv");


Nit = dim(data)[1];

I0 = c(1:Nit);

I1 = c(1:Nit);

I2 = c(1:Nit);

pH = 1.5574*10^(-4)

x = c(1:Nit);
cx = c(1:Nit);

a = c(1:Nit);

b = c(1:Nit);

for(i in 1:Nit)
{
  print(c(i,data$Peptide[i]))
  NEH = data$NEH[i];
  
  I0_0 = data$I0.0[i];#0.486729062;
  
  I1_0 = data$I1.0[i];#0.342671807;
  
  I2_0 = data$I2.0[i];#0.133438788;
  
  x[i] = 0.046;#0.0001 + i * 0.35 / Nit;
  
  I0[i] = data$I0.11[i];# 0.148189925 #I0_t(I0_0, x[i], NEH);
  
  I1[i] = data$I1.11[i];#0.238227398 #I1_t(I0_0, I1_0, NEH, x[i]);
  
  I2[i] = data$I2.11[i];# 0.264838501 #I2_t(I0_0, I1_0, I2_0, x[i], NEH);
  
  
  delta_A2_t = I2[i] / I0[i] - I2_0 / I0_0;
  
  delta_A1_t = I1[i] / I0[i] - I1_0 / I0_0;
  
  #print(c(i, delta_A1_t, delta_A2_t));
  
  A = delta_A2_t - (I1_0 / I0_0) * pH * delta_A1_t - (I1[i] / I0[i]) * (1 - pH) * delta_A1_t;
  
  A = A + (-delta_A1_t^2 * (2*pH - 1) + delta_A1_t * ( 2 * pH - 1) / (1 - pH) ) / 2;
  
  B = -delta_A2_t * (1 - pH) + (I1_0 / I0_0) * pH * delta_A1_t * 2 * ( 1 - pH);
  
  B = B + (1 - 2 * pH) * I1[i] / I0[i] * ( 1 - pH) * delta_A1_t;
  
  B = B + ( delta_A1_t^2 * (1 - pH) * ( 2* pH - 1) + delta_A1_t^2 * 2 * pH * (1 - pH) - delta_A1_t * 2 * pH ) / 2;
  
  # format(print(c(I0[i], I1[i], I2[i], x[i], A, B, -B / A,'test')), digits=4);
  
  a[i] = A; b[i] = B;
  cx[i]=-B / A
}

df<- data.frame(data$Peptide, a,b,cx,x)
write.csv(df, "temp_res2.csv", row.names = FALSE)
