# Define your signal function (replace with your specific function)
signal_func <- function(x) {
  sin(2 * pi * x)
}

# Set desired SNR and number of data points
snr <- 20 # In dB
n <- 100

# Generate data points
x <- seq(0, 1, length = n)

# Calculate noise standard deviation based on SNR
sigma <- sqrt(var(signal_func(x)) / (10^(snr / 10)))

# Generate Gaussian noise
noise <- rnorm(n, mean = 0, sd = sigma)

# Add noise to signal and create final data
data <- signal_func(x) + noise

# Plot the data
plot(x, data, main = "Simulated Data with SNR")