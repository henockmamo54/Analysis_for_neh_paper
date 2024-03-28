


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




# Set desired SNR and number of data points
snr <- 80 # In dB
n <- 110

# Generate data points 
# x=c(0.4581,0.3416, 0.1434,0.125,0.06,0.03)
data=read.csv("C:\\Workplace\\C#\\Test\\Test_test\\ConsoleApp2\\bin\\Debug\\0_01.csv");
index=10;
x=c(data$I0[index],data$I1[index],data$I2[index],
    data$I3[index],data$I4[index],data$I5[index]);

x=x/sum(x)

# Calculate noise standard deviation based on SNR
sigma <- sqrt(var((x)) / (10^(snr / 10)))
print(c('sigma',sigma))

# Generate Gaussian noise
# noise <- rnorm(length(x), mean = 0, sd = sigma)
temp_noise=generateNoisyData(x,snr = snr)
noise=temp_noise[1:6]
sigma=temp_noise[7]
c_snr=temp_noise[8]


# Add noise to signal and create final data
# data <- signal_func(x) + noise
y <- (x) + noise

# Plot the data
plot(x, y, main = "Simulated Data with SNR")
lines(x, x,col ='red', main = "Simulated Data with SNR")

computed_snr <- 10*log10(var(x) / var(noise)) #sqrt(sum(x^2) / sum(noise^2))#
print(paste("Specified SNR:", snr))
print(paste("Computed SNR:", computed_snr))

y=y/sum(y)
print(data$NEH[index])
Compute_neh_pxt(data$M0[index],data$M1[index],data$M2[index], 
                y[1],y[2],y[3]);

