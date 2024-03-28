


NEH_Compute <- function(A0, A1, A0_0, A1_0, pX)
{
  pH = 1.5574*10^(-4)
  NEH = (1 - pH - pX)*(1 - pH) / pX;
  NEH = NEH * (A1 / A0 - A1_0 / A0_0);
  return(NEH)
}

generateNoisyData<- function(x,snr){
  sigma <- sqrt(var((x)) / (10^(snr / 10)))
  noise <- rnorm(length(x), mean = 0, sd = sigma)
  data <- (x) + noise
  computed_snr <- 10*log10(var(x) / var(noise))
  while (abs(snr-computed_snr)>0.01) {
    sigma <- sqrt(  (var((x))+ mean(x)*mean(x)) / (10^(snr / 10)))
    noise <- rnorm(length(x), mean = 0, sd = sigma)
    computed_snr <- 10*log10(var(x) / var(noise))
  }
  
  
  data <- (x) + noise
  return(c(noise,sigma,computed_snr))
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
 

# Set desired SNR and number of data points
# snr <- 120 # In dB 

# Generate data points 
# x=c(0.4581,0.3416, 0.1434,0.125,0.06,0.03)
data=read.csv("C:\\Workplace\\C#\\Test\\Test_test\\ConsoleApp2\\bin\\Debug\\0_05.csv");
# data=read.csv("C:\\Workplace\\C#\\Test\\Test_test\\ConsoleApp2\\bin\\Debug\\0_03.csv");
# data=read.csv("C:\\Workplace\\C#\\Test\\Test_test\\ConsoleApp2\\bin\\Debug\\0_01.csv");
pw=0.05
# pw=0.03
# pw=0.01

set.seed(40)  
neh_dif=list()
pxt_dif=list()
for (snr in c(120,90,60,40,30,10)){
  print(c("==> snr",snr))
  for (index in seq(1:100)){
    # print(c("index ",index))
    x=c(data$I0[index],data$I1[index],data$I2[index],
        data$I3[index],data$I4[index],data$I5[index]);
    x=x/sum(x)
    
    # Calculate noise standard deviation based on SNR
    sigma <- sqrt(var((x)) / (10^(snr / 10)))
    temp_neh_dif=list()
    temp_pxt_dif=list()
    for(sim_index in seq(1,10)){
      # Generate Gaussian noise
      # noise <- rnorm(length(x), mean = 0, sd = sigma)
      temp_noise=generateNoisyData(x,snr = snr)
      noise=temp_noise[1:6]
      sigma=temp_noise[7]
      c_snr=temp_noise[8]
      
      
      # Add noise to signal and create final data
      # data <- signal_func(x) + noise
      y <- (x) + noise
      
      
      y=y/sum(y) 
      res=Compute_neh_pxt(data$M0[index],data$M1[index],data$M2[index], 
                      y[1],y[2],y[3]);
      temp_neh_dif=append(temp_neh_dif,100*abs(data$NEH[index]- res[2])/data$NEH[index])
      temp_pxt_dif=append(temp_pxt_dif,100*abs(pw-res[1])/pw)
    }
    
    # print(c( abs(pw-res[1])/pw, 
    #          abs(data$NEH[index]- res[2])/data$NEH[index])
    #       )
    neh_dif=append(neh_dif, median(unlist(temp_neh_dif)) )
    pxt_dif=append(pxt_dif, median(unlist(temp_pxt_dif)) )
  }
  
  print(c(median(unlist(neh_dif)),median(unlist(pxt_dif))))
}
# 
df <- data.frame(unlist(neh_dif),unlist(pxt_dif))
write.csv(df, "sim_pw_05.csv", row.names = FALSE)
# write.csv(df, "sim_pw_03.csv", row.names = FALSE)
# write.csv(df, "sim_pw_01.csv", row.names = FALSE)
print('done!')

