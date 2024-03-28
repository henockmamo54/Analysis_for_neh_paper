
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
 

snr=90;
# pw=0.05;
pw=0.03;
# pw=0.01;
# data=read.csv("C:\\Workplace\\C#\\Test\\Test_test\\ConsoleApp2\\bin\\Debug\\0_05.csv");
data=read.csv("C:\\Workplace\\C#\\Test\\Test_test\\ConsoleApp2\\bin\\Debug\\0_03.csv");
# data=read.csv("C:\\Workplace\\C#\\Test\\Test_test\\ConsoleApp2\\bin\\Debug\\0_01.csv");

file_neh_dif=list()
file_pxt_dif=list()

for(n_snr in c(90,50,30,10)){
  print(n_snr)
  snr_neh_dif=list()
  snr_pxt_dif=list()
  for(index in seq(1,100)){
    print(c('\t\t=>',index))
    temp_neh_dif=list()
    temp_pxt_dif=list()
    for(sim_index in seq(1,10)){
      x=c(data$I0[index],data$I1[index],data$I2[index],
          data$I3[index],data$I4[index],data$I5[index]);
      
      temp_noise=generateNoisyData(x,snr = snr)
      noise=temp_noise[1:6]
      sigma=temp_noise[7]
      c_snr=temp_noise[8]
      
      
      # #add noise and normalize
      x=x+noise
      sum_x_noise=sum(x)
      x=x/sum_x_noise
      
      
      # Compute_neh_pxt(data$M0[index],data$M1[index],data$M2[index], 
      #                 data$I0[index],data$I1[index],data$I2[index]);
      
      
      res=Compute_neh_pxt(data$M0[index],data$M1[index],data$M2[index], 
                          x[1],x[2],x[3]);
      
      # print(c(data$NEH[index],res[2],c_snr,snr))
      temp_neh_dif=append(temp_neh_dif,abs(data$NEH[index]-res[2])/data$NEH[index])
      temp_pxt_dif=append(temp_pxt_dif,abs(pw-res[1])/pw)
    }
    
    # print(c(median(unlist(temp_pxt_dif)),median(unlist(temp_neh_dif))))
    snr_neh_dif=append(snr_neh_dif, median(unlist(temp_neh_dif)) )
    snr_pxt_dif=append(snr_pxt_dif, median(unlist(temp_pxt_dif)) )
  }
  
  file_neh_dif=append(file_neh_dif,snr_neh_dif)
  file_pxt_dif=append(file_pxt_dif,snr_pxt_dif)
}


df <- data.frame(unlist(file_neh_dif),unlist(file_pxt_dif))
# write.csv(df, "sim_pw_05.csv", row.names = FALSE)
write.csv(df, "sim_pw_03.csv", row.names = FALSE)
# write.csv(df, "sim_pw_01.csv", row.names = FALSE)
print('done!')



