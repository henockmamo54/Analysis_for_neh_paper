
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
  
  return (c(new_pxt, new_neh));
  
}


## sample
#==========================

NEH = 23;
x=0.046;
# theoretical unlabeled mass isotopomers
I0_0 = 0.4581;
I1_0 = 0.3416;
I2_0 = 0.1434; 
# labeled mass isotopomers
I0 =  I0_t(I0_0, x, NEH);
I1 = I1_t(I0_0, I1_0, NEH, x);
I2 = I2_t(I0_0, I1_0, I2_0, x, NEH);

vals=Compute_neh_pxt(I0_0,I1_0,I2_0,I0,I1,I2); # computes the px(t) value and 
                                               #the corresponding NEH value

print(c("px(t)",vals[1],"NEH",vals[2]))




