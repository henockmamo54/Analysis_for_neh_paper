
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


Profile_Peptide4 <- function() # theoretical
{
  
  data=read.csv("C:\\Workplace\\Python\\AnalysisForThePaper\\NEH\\second_round\\prepare_Data\\data.csv");
  Nit = dim(data)[1];
  pH = 1.5574*10^(-4) 
  
  I0 = c(1: (Nit*100));
  I1 = c(1:(Nit*100));
  I2 = c(1:(Nit*100));  
  
  x = c(1:(Nit*100));
  cx = c(1:(Nit*100));
  
  a = c(1:(Nit*100));  
  b = c(1:(Nit*100));
  
  c_neh=c(1:(Nit*100));
  t_neh=c(1:(Nit*100));
  
  # pep=c(1:(Nit*100));
  # I0_0 = c(1: (Nit*100));
  # I1_0 = c(1:(Nit*100));
  # I2_0 = c(1:(Nit*100)); 
  pep = c(rep("-1", length((Nit*100))));
  
  delta_1=c(1:(Nit*100));
  delta_2=c(1:(Nit*100));
  
  a1=c(1:(Nit*100));
  a2=c(1:(Nit*100));
  
  a_20=c(1:(Nit*100));
  a_21=c(1:(Nit*100));
  a_22=c(1:(Nit*100));
  a_23=c(1:(Nit*100));
  
  b1=c(1:(Nit*100));
  b2=c(1:(Nit*100));
  b3=c(1:(Nit*100));
  
  
  b_20=c(1:(Nit*100));
  b_21=c(1:(Nit*100));
  b_23=c(1:(Nit*100));
  b_22=c(1:(Nit*100));
  b_24=c(1:(Nit*100));
  b_25=c(1:(Nit*100)); 
  
  rmse=c(1:(Nit*100));
  
  
  
  for(i in 1:Nit)
  {
    print(i)
    for (j in 1:100)
    {
      
      x[i*100+j] =0.046#0.0001 + j * 0.1 / 100;# 0.046; #0.0001 + i * 0.35 / Nit;
      
      NEH = data$NEH[i];
      t_neh[i*100+j]=NEH;
      pep[i*100+j]=data$Peptide[i];
      # print(c("=>>>",NEH,"****",x[i+j]));
      
      I0_0 = data$M0[i]/100;# data$I0.0[i];
      
      I1_0 = data$M1[i]/100;#data$I1.0[i];
      
      I2_0 = data$M2[i]/100;#data$I2.0[i];
      
      # temp_sum= I0_0+I1_0+I2_0;
      # I0_0=I0_0/temp_sum;
      # I1_0=I1_0/temp_sum;
      # I2_0=I2_0/temp_sum;
      # 
      #================================
      
      # pep[i*100+j]=data$Peptide;
      # I0_0[i*100+j]=I0_0;
      # I1_0[i*100+j]=I1_0;
      # I2_0[i*100+j]=I2_0;
      #================================
      
      
      
      I0[i*100+j] =  I0_t(I0_0, x[i*100+j], NEH); #data$I0.11[i];#  I0_t(I0_0, x[i*100+j], NEH) +rnorm(1,0,0.005) ; # 
      
      I1[i*100+j] =  I1_t(I0_0, I1_0, NEH, x[i*100+j]); #data$I1.11[i];#I1_t(I0_0, I1_0, NEH, x[i*100+j])+rnorm(1,0,0.005) ; #
      
      I2[i*100+j] =  I2_t(I0_0, I1_0, I2_0, x[i*100+j], NEH); #data$I2.11[i];#I2_t(I0_0, I1_0, I2_0, x[i*100+j], NEH)+rnorm(1,0,0.005); # 
      
      # temp_sum=I0[i*100+j]+I1[i*100+j]+I2[i*100+j]
      # I0[i*100+j]=I0[i*100+j]/temp_sum;
      # I1[i*100+j]=I1[i*100+j]/temp_sum;
      # I2[i*100+j]=I2[i*100+j]/temp_sum;
        
      
      
      delta_A2_t = I2[i*100+j] / I0[i*100+j] - I2_0 / I0_0;
      
      delta_A1_t = I1[i*100+j] / I0[i*100+j] - I1_0 / I0_0;
      
      delta_2[i*100+j]=delta_A2_t;
      delta_1[i*100+j]=delta_A1_t;
      
      #print(c(i, delta_A1_t, delta_A2_t));
      
      A = delta_A2_t - (I1_0 / I0_0) * pH * delta_A1_t - (I1[i*100+j] / I0[i*100+j]) * (1 - pH) * delta_A1_t;
      a1[i*100+j]=A
      
      A = A + (-delta_A1_t^2 * (2*pH - 1) + delta_A1_t * ( 2 * pH - 1) / (1 - pH) ) / 2;
      a2[i*100+j]=(-delta_A1_t^2 * (2*pH - 1) + delta_A1_t * ( 2 * pH - 1) / (1 - pH) ) / 2
      
      a_20[i*100+j]=delta_A2_t;
      a_21[i*100+j]=-(I1[i*100+j] / I0[i*100+j]) * delta_A1_t;
      a_22[i*100+j]=0.5*delta_A1_t^2;
      a_23[i*100+j]=-0.5*delta_A1_t;
      
      #=============================================B==================================
      B = -delta_A2_t * (1 - pH) + (I1_0 / I0_0) * pH * delta_A1_t * 2 * ( 1 - pH);
      b1[i*100+j]=B
      
      B = B + (1 - 2 * pH) * I1[i*100+j] / I0[i*100+j] * ( 1 - pH) * delta_A1_t;
      b2[i*100+j]=(1 - 2 * pH) * I1[i*100+j] / I0[i*100+j] * ( 1 - pH) * delta_A1_t;
      
      B = B + ( delta_A1_t^2 * (1 - pH) * ( 2* pH - 1) + delta_A1_t^2 * 2 * pH * (1 - pH) - delta_A1_t * 2 * pH ) / 2;
      b3[i*100+j]=( delta_A1_t^2 * (1 - pH) * ( 2* pH - 1) + delta_A1_t^2 * 2 * pH * (1 - pH) - delta_A1_t * 2 * pH ) / 2;
      # format(print(c(I0[i*100+j], I1[i*100+j], I2[i*100+j], x[i*100+j], A, B, -B / A, NEH_Compute(I0[i*100+j],I1[i*100+j], I0_0,I1_0,pX = x[i*100+j] ),NEH)), digits=4);
      
      b_20[i*100+j]= -delta_A2_t * (1 - pH); #<-
      b_21[i*100+j]=(I1_0 / I0_0) * pH * delta_A1_t * 2 * ( 1 - pH);
      b_22[i*100+j]=(1 - 2 * pH) * I1[i*100+j] / I0[i*100+j] * ( 1 - pH) * delta_A1_t;#<-
      b_23[i*100+j]=(delta_A1_t^2 * (1 - pH) * ( 2* pH - 1)/2); #<-
      b_24[i*100+j]=(delta_A1_t^2 * 2 * pH * (1 - pH))/2;
      b_25[i*100+j]=( delta_A1_t * 2 * pH )/2 ;
      
      #=====================================delete=====================
      
      
      # delta_A2_t = I2_t(I0_0, I1_0, I2_0, x[i*100+j], NEH) / I0_t(I0_0, x[i*100+j], NEH) - I2_0 / I0_0;
      # 
      # delta_A1_t = I1_t(I0_0, I1_0, NEH, x[i*100+j]) / I0_t(I0_0, x[i*100+j], NEH) - I1_0 / I0_0;
      #=======================================end delete ===================
      c= delta_A2_t - (I1[i*100+j] / I0[i*100+j])*delta_A1_t + (delta_A1_t^2)*0.5
      temp=1 + (0.5/ ( (delta_A2_t/delta_A1_t) - (I1[i*100+j] / I0[i*100+j]) + delta_A1_t/2 - 0.5))
      temp2=1 + (0.5/ ( (delta_A2_t/delta_A1_t) - (I1[i*100+j] / I0[i*100+j]) + delta_A1_t/2 - 0.5))
      
      temp2=1 + (0.5/ ( (delta_A2_t/delta_A1_t) - delta_A1_t/2 - I1_0 / I0_0 -  0.5))
      
      #temp2= 1 + (0.5/ ( 0.8 - delta_A1_t/2 - I1_0 / I0_0 -  0.5))
      
      #temp2= 1+((0.5*delta_A1_t)/ (delta_A2_t- (delta_A1_t/2 ) + (delta_A1_t^2)/2 ));
      
      # temp=1 + (0.5/ ( (delta_A2_t/delta_A1_t) - (I1_t(I0_0, I1_0, NEH, x[i*100+j]) / I0_t(I0_0, x[i*100+j], NEH)) + delta_A1_t/2 - 0.5))
      # temp=1 + (0.5/ ( (delta_A2_t/delta_A1_t) - (I1_t(I0_0, I1_0, NEH, x[i*100+j]) / I0_t(I0_0, x[i*100+j], NEH)) + delta_A1_t/2 - 0.5))
      
      # temp=1 + (0.5/ ( (delta_A2_t/delta_A1_t) - (  I1_t(I0_0, I1_0, NEH, x[i*100+j])  / I0_t(I0_0, x[i*100+j], NEH)) + delta_A1_t/2 - 0.5))
      # temp=c/(c -delta_A1_t/2 )
      # temp= 1 + ((delta_A1_t/2)/(c-delta_A1_t/2))
      print(c(x[i*100+j],-B/A,temp,temp2))
      # print(c('**',-B,c))
      # print(c('**',A,c - delta_A1_t/2))
      # print(c('**',-B/A,c/(c - delta_A1_t/2)))
      #=============================================B==================================
      
      a[i*100+j] = A; 
      b[i*100+j] = B;
      c_neh[i*100+j]=NEH_Compute(I0[i*100+j],I1[i*100+j], I0_0,I1_0,pX = x[i*100+j] );
      cx[i*100+j]=-B / A;
      
      rmse[i*100+j]= sqrt( 
        ((I0[i*100+j] - I0_t(I0_0, cx[i*100+j], NEH) )^2 +
           (I1[i*100+j] - I1_t(I0_0, I1_0, NEH, cx[i*100+j])  )^2 +
           (I2[i*100+j] - I2_t(I0_0, I1_0, I2_0, cx[i*100+j], NEH) )^2 
        )/3);
      
      
    }
    
    # if (i==5){break;}
    
  }
  
  
  print(c( length(a),length(b),length(cx),length(x),length(t_neh),length(c_neh)   ))
  
  df <- data.frame(pep[101:length(a)],
                   # pep[101:length(pep)] ,
                   # I0_0[101:length(I0_0)] ,
                   # I1_0[101:length(I1_0)] ,
                   # I2_0[101:length(I2_0)] ,
                   a[101:length(a)],
                   b[101:length(b)],
                   cx[101:length(cx)],
                   x[101:length(x)],
                   t_neh[101:length(t_neh)],
                   c_neh[101:length(c_neh)],
                   b1[101:length(a)],b2[101:length(a)],b3[101:length(a)],
                   a1[101:length(a)],a2[101:length(a)],
                   delta_1[101:length(a)],delta_2[101:length(a)],
                   rmse[101:length(a)],
                   a_20[101:length(a)],a_21[101:length(a)],a_22[101:length(a)],a_23[101:length(a)],
                   b_20[101:length(a)],b_21[101:length(a)],b_22[101:length(a)],b_23[101:length(a)],b_24[101:length(a)],b_25[101:length(a)],
                   I0[101:length(a)], I1[101:length(a)]
  )
  
  write.csv(df, "res_theo.csv", row.names = FALSE)
  print('done!')
  
}

