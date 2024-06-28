

remove(list = ls())


NEH_Compute <- function(A0, A1, A0_0, A1_0, pX)
{
  pH = 1.5574 * 10 ^ (-4)
  NEH = (1 - pH - pX) * (1 - pH) / pX
  
  NEH = NEH * (A1 / A0 - A1_0 / A0_0)
  
  return(NEH)
}


generateNoisyData <- function(x, snr) {
  sigma <- sqrt((var((x)) + mean(x) * mean(x)) / (10 ^ (snr / 10)))
  noise <- rnorm(length(x), mean = 0, sd = sigma)
  data <- (x) + noise
  computed_snr <-
    10 * log10((var((x)) + mean(x) * mean(x)) / var(noise))
  while (abs(snr - computed_snr) > 0.01) {
    sigma <- sqrt((var((x)) + mean(x) * mean(x)) / (10 ^ (snr / 10)))
    noise <- rnorm(length(x), mean = 0, sd = sigma)
    computed_snr <-
      10 * log10((var((x)) + mean(x) * mean(x)) / var(noise))
  }
  
  
  
  
  data <- (x) + noise
  return(c(noise, sigma, computed_snr))
}

Compute_neh_pxt <- function(I0_0, I1_0, I2_0, I0, I1, I2) {
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


# Set desired SNR and number of data points
# snr <- 120 # In dB


test_path = "C:\\Workplace\\Python\\AnalysisForThePaper\\NEH\\Final\\C# simulation data gen\\ConsoleApp2\\bin\\Debug\\"
data = read.csv(paste(test_path, "0_035.csv", sep = ''))


pw = 0.035
snr_levels = seq(1, 40, 2)

# set.seed(40)
set.seed(60)
neh_dif = list()
pxt_dif = list()
rss = list()
rss2 = list()



final_res_colnames <-
  c(    "Protein",    "Peptide" ,    "Charge" ,    "NEH"   ,    "M0"   ,
    "M1"   ,    "M2" ,    "M3"  ,    "M4"   ,    "M5"   ,    "M0_1"   ,
    "M1_1"   ,    "M2_1" ,    "M3_1"  ,    "M4_1"   ,    "M5_1"   ,
    "I0"     ,    "I1"   ,    "I2"  ,    "I3"   ,    "I4"   ,    "I5"   ,
    "pxt",    "SNR"
  )
final_res <-
  data.frame(matrix(ncol = length(final_res_colnames), nrow = 0))
colnames(final_res) <-
  final_res_colnames #generate dataframe to store parameters and calculate values


# for (snr in c(120,90,60,40,30,10)){
for (snr in snr_levels) {
  print(c("==> snr", snr))
  for (index in seq(1:100)) {
    print(index)
    # print(c("index ",index))
    x = c(data$I0[index], data$I1[index], data$I2[index],
          data$I3[index], data$I4[index], data$I5[index])
    
    x = x / sum(x)
    
    
    # Calculate noise standard deviation based on SNR
    # sigma <- sqrt(var((x)) / (10^(snr / 10)))
    temp_neh_dif = list()
    temp_pxt_dif = list()
    temp_rss = list()
    temp_rss2 = list()
    for (sim_index in seq(1, 10)) {
      # Generate Gaussian noise
      # noise <- rnorm(length(x), mean = 0, sd = sigma)
      temp_noise = generateNoisyData(x, snr = snr)
      noise = temp_noise[1:6]
      sigma = temp_noise[7]
      c_snr = temp_noise[8]
      
      
      # Add noise to signal and create final data
      # data <- signal_func(x) + noise
      y <- (x) + noise
      
      
      final_res[nrow(final_res) + 1, ] = list(
        data$Protein[index],
        data$Peptide[index],
        data$Charge[index],
        data$NEH[index],
        data$M0[index] * 100,
        data$M1[index] * 100,
        data$M2[index] * 100,
        data$M3[index] * 100,
        data$M4[index] * 100,
        data$M5[index] * 100,
        data$M0[index] * 100,
        data$M1[index] * 100,
        data$M2[index] * 100,
        data$M3[index] * 100,
        data$M4[index] * 100,
        data$M5[index] * 100,
        y[1],
        y[2],
        y[3],
        y[4],
        y[5],
        y[6],
        0.035,
        snr
      )
      
    }
    
  }
  
}

write.csv(final_res, "0_035_with_noise.csv", row.names = FALSE)
print('done!')
