
NEH_Compute <- function(A1, A2, A1_0, A2_0, pX)
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


Profile_Peptide <- function()
{

   NEH = 23;

   Nit = 100;

   I0_0 = 0.4581;

   I1_0 = 0.3416;

   I2_0 = 0.1434;
   
   I0 = c(1:Nit);

   I1 = c(1:Nit);

   I2 = c(1:Nit);

   pH = 1.5574*10^(-4)

   x = c(1:Nit);

   a = c(1:Nit);

   b = c(1:Nit);

   for(i in 1:Nit)
   {
      x[i] = 0.0001 + i * 0.35 / Nit;

      I0[i] =  I0_t(I0_0, x[i], NEH);

      I1[i] = I1_t(I0_0, I1_0, NEH, x[i]);

      I2[i] = I2_t(I0_0, I1_0, I2_0, x[i], NEH);


      delta_A2_t = I2[i] / I0[i] - I2_0 / I0_0;

      delta_A1_t = I1[i] / I0[i] - I1_0 / I0_0;

      #print(c(i, delta_A1_t, delta_A2_t));

      A = delta_A2_t - (I1_0 / I0_0) * pH * delta_A1_t - (I1[i] / I0[i]) * (1 - pH) * delta_A1_t;

      A = A + (-delta_A1_t^2 * (2*pH - 1) + delta_A1_t * ( 2 * pH - 1) / (1 - pH) ) / 2;

      B = -delta_A2_t * (1 - pH) + (I1_0 / I0_0) * pH * delta_A1_t * 2 * ( 1 - pH);

      B = B + (1 - 2 * pH) * I1[i] / I0[i] * ( 1 - pH) * delta_A1_t;

      B = B + ( delta_A1_t^2 * (1 - pH) * ( 2* pH - 1) + delta_A1_t^2 * 2 * pH * (1 - pH) - delta_A1_t * 2 * pH ) / 2;

      format(print(c(I0[i], I1[i], I2[i], x[i], A, B, -B / A)), digits=4);

      a[i] = A; b[i] = B;
   }
   

   plot(x, I0, type="l", lwd=2, xlab = "Deuterium Enrichment", ylab="Relative Abundance");

   lines(x, I1, col="blue", lwd=2);

   lines(x, I2, col="purple", lwd=2);

   legend("topright",
   c(expression(paste("I"[0], "(t)")), expression(paste("I"[1], "(t)")),
   expression(paste("I"[2], "(t)"))), col=c("black", "blue", "purple"), lwd=c(2,2,2));

   legend("top", bty="n", "IQDAGLVLADALR");

   plot(a, b);
}


Compute_PXT <- function(data_peptides)
{
   pH = 1.5574*10^(-4);

   Npeptides = length(data_peptides$Peptide);

   a = c(1:Npeptides);

   b = c(1:Npeptides);

   PXT = c(1:Npeptides);

   Nprint = 5;

   #Nprint = Npeptides;

   #for(i in 1:Npeptides)
   for(i in 1:Nprint)
   {

      I0_0 = data_peptides$I0.0[i];     I1_0 = data_peptides$I1.0[i];    I2_0 = data_peptides$I2.0[i];

      I0_0 = data_peptides$M0[i];       I1_0 = data_peptides$M1[i];      I2_0 = data_peptides$M2[i];

      I0_t = data_peptides$I0.11[i];    I1_t = data_peptides$I1.11[i];   I2_t = data_peptides$I2.11[i];


      print(noquote(c(I0_0, I1_0, I2_0, I0_t, I1_t, I2_t)), digits=4 );


      #I0_0 = 0.4581;         I1_0 = 0.3416;           I2_0 = 0.1434;

      #I0_t = 0.330441034;    I1_t = 0.355134877;      I2_t = 0.201609870;

      delta_A2_t = I2_t / I0_t - I2_0 / I0_0;

      delta_A1_t = I1_t / I0_t - I1_0 / I0_0;

      A = delta_A2_t - (I1_0 / I0_0) * pH * delta_A1_t - (I1_t / I0_t) * (1 - pH) * delta_A1_t;

      print(c("A first : ", A));

      A = A + (-delta_A1_t^2 * (2*pH - 1) + delta_A1_t * ( 2 * pH - 1) / (1 - pH) ) / 2;
      
      print(c("Second : ", (-delta_A1_t^2 * (2*pH - 1) + delta_A1_t * ( 2 * pH - 1) / (1 - pH) ) / 2));


      B = -delta_A2_t * (1 - pH) + (I1_0 / I0_0) * pH * delta_A1_t * 2 * ( 1 - pH);

      B = B + (1 - 2 * pH) * I1_t / I0_t * ( 1 - pH) * delta_A1_t;

      B = B + ( delta_A1_t^2 * (1 - pH) * ( 2* pH - 1) + delta_A1_t^2 * 2 * pH * (1 - pH) - delta_A1_t * 2 * pH ) / 2;


      a[i] = A; b[i] = B; PXT[i] = - B / A;

      print(noquote(format(c(i, A, B, -B / A), digits=4 )));

      plot(density(PXT));

   }


   plot(density(PXT[1:Nprint]), xlim=c(-20, 20));

   print(median(PXT[1:Nprint]));

   plot(a[1:Nprint], b[1:Nprint]);

   print(summary(lm(b[1:Nprint] ~  a[1:Nprint] + 0)));


}

Simulate_Rate <- function(k_deg = 0.27, pW = 0.03, SD = 0.01)
{

   NEH = 23;

   Nit = 100;

   I0_0 = 0.4581;

   I1_0 = 0.3416;

   I2_0 = 0.1434;

   pH = 1.5574*10^(-4);

   t1 <- c(0, 3, 5, 7, 14, 21);

   I0_asymp = I0_0 * (1 - pW/(1 - pH))^NEH;

#   print(c(I0_0, pW, pH, NEH));

   I0_t = I0_asymp + (I0_0 - I0_asymp) * exp(-k_deg * t1);

   temp = (1/NEH) * log((1 - pW/(1-pH))^NEH + (1 - (1 - pW/(1-pH))^NEH) * exp(-k_deg * t1));

   pX_t = (1 - exp(temp)) * (1 - pH);

   pX_t = 1 - ((1 - pW/(1-pH))^NEH + (1 - (1 - pW/(1-pH))^NEH) * exp(-k_deg * t1))^(1 / NEH) ;

   pX_t = pX_t * (1 - pH);

   I0_X1 = I0_0 * (1 - pX_t/(1 - pH)) ^NEH;

#   print(I0_X1 - I0_t);

#   print(I0_asymp);

   print(I0_t);

   print(pX_t);

 I00_t = I0_t;

 for(i in 1:10)
{
   noise = rnorm(length(t1), 0.0, SD);

   I0_t = I00_t + noise;

   mod1 <- nls(I0_t ~ I0_asymp + (I0_0 - I0_asymp) * exp(-k * t1),
           start = list(k = 0.001));

   plot(t1, I0_t, type="l", lwd=2);

   plot(t1, pX_t, type="l", lwd=2);

   k_nls = coef(mod1)[[1]];

 #  I0_X1 = I0_X1 + noise;

#   mod2 <- nls(I0_X1 ~ I0_asymp + (I0_0 - I0_asymp) * exp(-k * t1),
 #          start = list(k = 0.001));
 

 #  print(c(k_nls, coef(mod2)[[1]]));


 #  print(k_nls - coef(mod2)[[1]]);



   #### Reconstract the I0(t) from A1(t) and A2(t) #####
   
   I1_t <- c(1:length(pX_t));

   I2_t <- c(1:length(pX_t));

   for(i in 1:length(pX_t))
   {
      I1_t[i] = I1_t(I0_0, I1_0, NEH, pX_t[i]);

      I2_t[i] = I2_t(I0_0, I1_0, I2_0, pX_t[i], NEH);
   }

   I1_t = I1_t + noise;

   beta = I0_t / (I0_t + I1_t);

   alpha = ( 1 / beta - 1 - I1_0 / I0_0) * (1 - pH) / NEH;

   pX_computed = (1 - pH) * alpha / (1 + alpha);

   I0_computed = I0_0 * (1 - pX_computed/(1-pH))^NEH;

   #print(I0_computed);

   #print(pX_t - pX_computed);

   mod3 = nls(I0_computed ~ I0_asymp + (I0_0 - I0_asymp) * exp(-k * t1),
           start = list(k = 0.001));

   #print(coef(mod3)[[1]]);

   print(c((coef(mod3)[[1]] - k_deg),  (coef(mod1)[[1]] - k_deg)));
 }


}

rlaplace <- function (n = 1, m = 0, s = 1)
{
    if (any(s <= 0)) 
        stop("s must be positive")
    q <- runif(n)
    ifelse(q < 0.5, s * log(2 * q) + m, -s * log(2 * (1 - q)) + 
        m)
}
          

#
#
#   A function to simulate rate constants
#
#    tm <- Rate_Constant_Simulation_A0A1 (NTrials = 1000);
#
#
Rate_Constant_Simulation_A0A1 <- function(NTrials = 2000, meanLog = -2.5988987, sdLog = 0.7316477,
   minRate = 0.01, maxRate = 0.6, noise_mean = -0.00375364,  noise_scale = exp(-4.140167273), BWE = 0.051)
{
   # do a number of iteration
   # with different starting values
   # of I_0, and I_asymp
   #

   tt = c(0, 3, 5, 7, 14, 21);

   pH = 0.00015574 / (0.99984426 + 0.00015574);

   nDays = length(tt);


   # Simulates the noise which will be
   # specific for each rate constant
   # and time point, but the same
   # for all timepoint sets

   Corrs_I0 <- c(1:NTrials);

   Corrs_I1 <- c(rep(0, NTrials));

   # first (0 days), and last (nEndLabeling) days are fixed

   N_Peptide_Length = seq(1:NTrials);

   I0_0 = seq(1:NTrials);

   I1_0 = seq(1:NTrials);

   I_asymp = seq(1:NTrials);

   K_Deg_Rates = seq(1:NTrials);

   K_I0_Deg = seq(1:NTrials);
   
   K_I1_Deg = seq(1:NTrials);

   Determ_I1 = seq(1:NTrials);

   Determ_I0 = seq(1:NTrials);

   N_Exch = seq(1:NTrials);

   y = seq(1:nDays);

   pxt_computed = seq(1:nDays);

   pxt_actual = seq(1:nDays);


   kk = 1;

  # set.seed(122);

  
  for(k in 1:NTrials)
  {

     while(TRUE)
     {
         temp = rlnorm(1, meanLog, sdLog);

         if(temp >= minRate & temp <= maxRate)
         {
             K_Deg_Rates[kk] = temp;
  
             kk = kk + 1;

             break;
         }

      }

     if(kk == NTrials + 1)
     {
         break
     }

  } # for(k in 1:NTrials)

  print("Finished Simulating Rates");

 # set.seed(122);

  # this block computes N_Peptide_Length, and related N_Exch, I0_0, I_asymp
  for(kk in 1:NTrials)
  {
     #
     # peptide length
     #
     N_Peptide_Length[kk] = as.integer(rgamma(1, 9.094925059, 0.604736252) + 0.5);

     #these parameters for Gamma distribution
     # may be changed based on the distribution
     # of peptides in a particular dataset

     phi_N_Exch = 0.03306681;

     a0 = -2.99750; 

     a1 = 2.13329;

     mu_N_Exch = a0 + a1 * N_Peptide_Length[kk];

     alpha_N_Exch = 1 / phi_N_Exch;

     beta_N_Exch = alpha_N_Exch / mu_N_Exch;

     N_Exch[kk] = as.integer(rgamma(1, alpha_N_Exch, beta_N_Exch) + 0.5);

     #
     #  Finished N_Exch
     #

     #
     #compute I0_0[kk]
     #

     phi_beta = 259.9;

     theta_beta = 0.9079 - 0.0857 * N_Peptide_Length[kk];

     mu_beta = 0.5 + atan(theta_beta)/pi;

     alpha_beta = phi_beta * mu_beta;

     beta_beta = alpha_beta * (1 - mu_beta) / mu_beta;
  
     I0_0[kk] = rbeta(1, alpha_beta, beta_beta);

     #compute I_asymp[kk]
     I_asymp[kk] = I0_0[kk] * (1 - BWE / (1 - pH))^N_Exch[kk];

  }


  for(kk in 1:NTrials)
  {
     temp = 1;

     while(temp >= 0.95 || temp < 0.7)
     {
         temp = runif(1, 0, 1);

	 I1_0[kk] = temp;

	 temp = temp + I0_0[kk];
     }
  }


 for(k in 1:NTrials)
 {

   dsmall = 1000.0;

  # N_Exch[k] = 23;

   #K_Deg_Rates[k] = 0.27;

   #I0_0[k] = 0.4581;

   #I1_0[k] = 0.3416;


   asymp_temp = I_asymp[k];

   zero_temp  = I0_0[k];

   #asymp_temp = I0_0[k] * (1 - BWE / (1 - pH))^N_Exch[k];

   #print(c(I0_0[k], BWE, pH, N_Exch[k]));

   # I0_asymp = I0_0 * (1 - pW/(1 - pH))^NEH;

   Zero_Asymp_Diff = zero_temp - asymp_temp;

   #print(asymp_temp);

   #print(c("k = ", k, NTrials), quote = FALSE):

    rss = 0.0;

    n_it_true = 0;

   # this value of y is determined from the kth value of kdeg, and kth values
   # of zero_temp, and asymp_temp;

   y = asymp_temp + Zero_Asymp_Diff * exp(-K_Deg_Rates[k] * tt);

#
#  actual pxt_actual:
#
#

  # print(y);

   pxt_actual = 1 - (y / I0_0[k]) ^ (1/N_Exch[k])

   pxt_actual = pxt_actual * (1 - pH);

  # print(pxt_actual);


#
#  determine actual I1_t from I1_0, I0_0, and pxt_actual
#
  I1_t = c(1:length(pxt_actual));
  

   for(i in 1:length(pxt_actual))
   {
      I1_t[i] = I1_t(I0_0[k], I1_0[k], N_Exch[k], pxt_actual[i]);
   }


#   y_tilde = y + rlaplace(nDays, noise_mean, noise_scale);

#   I1_t = I1_t + rlaplace(nDays, noise_mean, noise_scale);

   #print(y);

   #print(I1_t);

   #noise_temp = rnorm(nDays, 0., 0.01);

   y_tilde = y + rnorm(nDays, 0., 0.01); #noise_temp;

   y_tilde = y + rlaplace(nDays, -0.00375364, exp(-4.140167273)); #noise_temp;

   #I1_t = I1_t + rnorm(nDays, 0., 0.01); #noise_temp;

#
#
#   Rate constant from y_tilde
#
#
   mod <- NULL;
   
   try( mod <- nls(y_tilde ~ asymp_temp + Zero_Asymp_Diff * exp(-a * tt), start = list(a = 0.1)), silent=TRUE);
   #mod <- nls(y ~ asymp_temp + Zero_Asymp_Diff * exp(-a * t), start = list(a = 0.1));
    #  print("Closed");

      if(!is.null(mod))
      {
          a = coef(mod)[[1]];

          #rss = rss + (a - K_Deg_Rates[k])^2 * dlnorm(K_Deg_Rates[k],meanLog, sdLog);
	  rss = rss + (a - K_Deg_Rates[k])^2 ;

          n_it_true = n_it_true + 1;

         #print(c(a, K_Deg_Rates[k]));

	 K_I0_Deg[k] = a;

         #   not clear about the last part
	 #Corrs_I0[k]= cor(y_tilde, asymp_temp + Zero_Asymp_Diff * exp(-a * tt)asymp_temp);

	 #Determ_I0[k] = sum((y_tilde - asymp_temp - Zero_Asymp_Diff * exp(-a * tt)asymp_temp)^2);

	 Determ_I0[k] = Determ_I0[k] / (length(y_tilde) - 1);

	 Determ_I0[k] = 1 - Determ_I0[k] / var(y_tilde);

        # print(a);


      }
      else
      {
           print(c("Bad inputs, i = ", k), quote = FALSE);
             #print(c(y_tilde, t, K_Deg_Rates[k]), quote = FALSE);

      }


#
#
#
#   Rate constant from I1_t and I0_t (which is y_tilde)
#
#   both I1_t and y_tilde have the noise added to them
#

    #beta = y_tilde / (y_tilde + I1_t);

    #alpha = ( 1 / beta - 1 - I1_0[k] / I0_0[k]) * (1 - pH) / N_Exch[k];

    # pX_computed = (1 - pH) * alpha / (1 + alpha);


    beta = I1_t / y - I1_0[k] / I0_0[k];

  #  beta = beta + rnorm(nDays, 0., 0.01);

    beta = beta + rlaplace(nDays, -0.00375364, exp(-4.140167273));

    alpha = (1 - pH) * beta;

    pX_computed = (1 - pH) * alpha / (N_Exch[k] + alpha);

    #print(pX_computed);

    I0_computed = I0_0[k] * (1 - pX_computed/(1-pH))^N_Exch[k];

    #print(c(asymp_temp, Zero_Asymp_Diff));

    #print(pX_computed);

    mod3 = nls(I0_computed ~ asymp_temp + Zero_Asymp_Diff  * exp(-k * tt),
           start = list(k = 0.001));


    K_I1_Deg[k] = coef(mod3)[[1]];

   # print(K_I1_Deg[k] - K_I0_Deg[k]);

    #plot(abs(K_I1_Deg[1:k] - K_Deg_Rates[1:k]), abs(K_Deg_Rates[1:k] - K_I0_Deg[1:k]),
    #xlim=c(0, 0.1), ylim=c(0., 0.1));

    #print(c((K_Deg_Rates[k] - K_I1_Deg[k]), (K_Deg_Rates[k] - K_I0_Deg[k])));
  }

   plot(K_I1_Deg, K_I0_Deg, xlab=expression(paste("Rate constants using A"[1], "(t)/", "A"[0], "(t)")),
       ylab=expression(paste("Rate constants using I "[0], "(t)")));

   results <- list();

   results$I1_0 = I1_0;
   results$K_Deg_Rates = K_Deg_Rates;
   results$I0_0 = I0_0;
   results$I_asymp = I_asymp;
   results$meanLog = meanLog;
   results$sdLog = sdLog;
   results$NTrials = NTrials;
   results$BWE = BWE;
   results$N_Peptide_Length = N_Peptide_Length;
   results$N_Exch = N_Exch;
   results$Corrs_I0 = Corrs_I0;
   results$Corrs_I1 = Corrs_I1;
   results$Determ_I0 = Determ_I0;
   results$Determ_I1 = Determ_I1;
   results$K_I1_Deg = K_I1_Deg;
   results$K_I0_Deg = K_I0_Deg;

  results
}