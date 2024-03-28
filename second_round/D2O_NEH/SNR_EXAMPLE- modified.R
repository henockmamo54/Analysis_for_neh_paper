
# x=c(0.4581,0.3416, 0.1434,0.5,0.6,0.18,0.75,0.95)
# 
# # Set desired SNR and number of data points
# snr <- 0 # In dB
# sigma <- sqrt(var(X) / (10^(snr / 10)))
# print(c("sigma",sigma));
# # Generate Gaussian noise
# noise <- rnorm(length(x), mean = 0, sd = sigma)
# 
# # Add noise to signal and create final data
# data <- x + noise
# 
# # Plot the data
# plot(x, data, main = "Simulated Data with SNR")
# lines(x, x,col ='red', main = "Simulated Data with SNR")
# 
# computed_snr <- sqrt(sum(x^2) / sum(noise^2))#10*log10(var(x) / var(noise))
# print(paste("Specified SNR:", snr))
# print(paste("Computed SNR:", computed_snr))


# Define your signal function (replace with your specific function)
# signal_func <- function(x) {
#   sin(2 * pi * x)
# }



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
  return(c(data,sigma,computed_snr))
}

# Set desired SNR and number of data points
snr <- -90 # In dB
n <- 110

# Generate data points
# x <- seq(0, 1, length = n)
# x=c(0.4581,0.3416, 0.1434,0.5,0.6,0.18,0.75,0.95)
x=c(0.4581,0.3416, 0.1434,0.125,0.06,0.03)
x=x/sum(x)

# Calculate noise standard deviation based on SNR
sigma <- sqrt(var((x)) / (10^(snr / 10)))
print(c('sigma',sigma))

# Generate Gaussian noise
noise <- rnorm(length(x), mean = 0, sd = sigma)

# Add noise to signal and create final data
# data <- signal_func(x) + noise
data <- (x) + noise

# Plot the data
plot(x, data, main = "Simulated Data with SNR")
lines(x, x,col ='red', main = "Simulated Data with SNR")

computed_snr <- 10*log10(var(x) / var(noise)) #sqrt(sum(x^2) / sum(noise^2))#
print(paste("Specified SNR:", snr))
print(paste("Computed SNR:", computed_snr))



