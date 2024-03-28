


NEH_Compute <- function(A0, A1, A0_0, A1_0, pX)
{
  #Computes the number of exchangeable hydrogen
  #A0 and A1 are normalized abundance of the ith mass isotopomer at labeling duration t,
  #A0_0 and A1_0 are the natural normalized abundance of the ith mass isotopomer
  #px - atomic percent excess in the deuterium
  pH = 1.5574 * 10 ^ (-4)
  NEH = (1 - pH - pX) * (1 - pH) / pX
  NEH = NEH * (A1 / A0 - A1_0 / A0_0)
  return(NEH)
}


I0_t <- function(I0_0, pX, NEH)
{
  # computes the normalized abundance of the mono isotope at labeling duration t,
  #I0_0 and I1_0 are the natural normalized abundance of the ith mass isotopomer
  #px - atomic percent excess in the deuterium
  #NEH - number of exchangeable hydrogen
  
  pH = 1.5574 * 10 ^ (-4)
  I0t = I0_0 * (1 - pX / (1 - pH)) ^ NEH
  return(I0t)
  
}

I1_t <- function(I0_0, I1_0, NEH, pX)
{
  # computes the normalized abundance of the first mass isotopomer at labeling duration t,
  #I0_0 and I1_0 are the natural normalized abundance of the ith mass isotopomer
  #px - atomic percent excess in the deuterium
  #NEH - number of exchangeable hydrogen
  pH = 1.5574 * 10 ^ (-4)
  I1t = (1 - pX / (1 - pH)) ^ (NEH - 1) * pX * NEH * I0_0 / (1 - pH) ^ 2
  I1t = I1t + I1_0 * (1 - pX / (1 - pH)) ^ NEH
  return(I1t)
  
}

I2_t <- function(I0_0, I1_0, I2_0, pX, NEH)
{
  # computes the normalized abundance of the second mass isotopomer at labeling duration t,
  #I0_0, I1_0 and I2_0 are the natural normalized abundance of the ith mass isotopomer
  #px - atomic percent excess in the deuterium
  #NEH - number of exchangeable hydrogen
  
  pH = 1.5574 * 10 ^ (-4)
  I0t = I0_t(I0_0, pX, NEH)
  I1t = I1_t(I0_0, I1_0, NEH, pX)
  I2t = I0t * (I2_0 / I0_0 - I1_0 * pH * NEH / (I0_0 * (1 - pH)))
  I2t = I2t + I0t * (pH / (1 - pH)) ^ 2 * NEH * (NEH + 1) / 2
  I2t = I2t - I0t * ((pX + pH) / (1 - pH - pX)) ^ 2 * NEH * (NEH + 1) / 2
  I2t = I2t + NEH * (pX + pH) * I1t / (1 - pH - pX)
  return(I2t)
  
}


Compute_neh_pxt <- function(I0_0, I1_0, I2_0, I0, I1, I2) {
  #Computes the number of exchangeable hydrogen and body water enrichment
  #I0_0 and I1_0 are the natural normalized abundance of the ith mass isotopomer
  #I0, I1 and I2 are normalized abundance of the ith mass isotopomer at labeling duration t,
  
  
  pH = 1.5574 * 10 ^ (-4)
  delta_A2_t = I2 / I0 - I2_0 / I0_0
  delta_A1_t = I1 / I0 - I1_0 / I0_0
  A = delta_A2_t - (I1_0 / I0_0) * pH * delta_A1_t - (I1 / I0) * (1 - pH) * delta_A1_t
  A = A + (-delta_A1_t ^ 2 * (2 * pH - 1) + delta_A1_t * (2 * pH - 1) / (1 - pH)) / 2
  B = -delta_A2_t * (1 - pH) + (I1_0 / I0_0) * pH * delta_A1_t * 2 * (1 - pH)
  B = B + (1 - 2 * pH) * I1 / I0 * (1 - pH) * delta_A1_t
  B = B + (delta_A1_t ^ 2 * (1 - pH) * (2 * pH - 1) + delta_A1_t ^ 2 * 2 * pH * (1 - pH) - delta_A1_t * 2 * pH) / 2
  new_pxt = -B / A
  new_neh = NEH_Compute(I0, I1, I0_0, I1_0, pX = new_pxt)
  
  return (c(new_pxt, new_neh))
  
  
}


## sample
#==========================

NEH = 23
x = 0.046

# theoretical unlabeled mass isotopomers
I0_0 = 0.4581
I1_0 = 0.3416
I2_0 = 0.1434

# labeled mass isotopomers
I0 =  I0_t(I0_0, x, NEH)
I1 = I1_t(I0_0, I1_0, NEH, x)
I2 = I2_t(I0_0, I1_0, I2_0, x, NEH)

# computes the px(t) value and
#the corresponding NEH value
vals = Compute_neh_pxt(I0_0, I1_0, I2_0, I0, I1, I2)

print(c("px(t)", vals[1], "NEH", vals[2]))


#===============================================
# simulation for multiple pep tides
#===============================================

data=read.csv("C:\\Workplace\\C++\\d2ome_restructure\\d2ome_restructure\\d2ome_GUI\\d2ome_GUI\\v2\\bin\\Debug\\Estimated_neh_numbers_utmb-liver_.csv")

pxt=c(rep(-1, dim(data)[1] ))
neh=c(rep(-1, dim(data)[1] ))

for (i in seq(1,(dim(data)[1]))){
  
  vals=Compute_neh_pxt(data$i0_0[i],data$i1_0[i],data$i2_0[i],
                       data$i0_31[i],data$i1_31[i],data$i2_31[i]); # computes the px(t) value and 
  #the corresponding NEH value
  
  print(c("px(t)",vals[1],"NEH",vals[2]))
  pxt[i]=vals[1];
  neh[i]=vals[2]
  
}

