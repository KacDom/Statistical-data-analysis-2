library(ggplot2)
library(dplyr)
library(gridExtra)
library(coda)

#Main function implementing Gibbs Sampler 
gibbs_sampler = function(num_of_samples){
    Rain = c(0)
    Cloudy = c(1)
    for (i in 2:num_of_samples){
        if (Rain[i-1] == 1){Cloudy[i] = sample(c(0,1), size = 1,  prob = c(0.556, 0.444))
        }else{Cloudy[i] = sample(c(0,1), size = 1, prob = c(0.952, 0.048))}
    
        if (Cloudy[i-1] == 1){Rain[i] = sample(c(0,1), size = 1, prob = c(0.185, 0.815))
        }else{Rain[i] = sample(c(0,1), size = 1, prob = c(0.784, 0.216))}
    }
    return(data.frame("Cloudy" = Cloudy, "Rain" = Rain))}


#Function for calculating relative frequency in every iteration of Gibbs sampling
count_relative_freq = function(gibbs_sample){
    relative_freq_Rain = c(0)
    relative_freq_Cloudy = c(0)
    for (i in 1:nrow(gibbs_sample)){
        relative_freq_Rain[i] = sum(gibbs_sample[1:i, 'Rain'])/i
        relative_freq_Cloudy[i] = sum(gibbs_sample[1:i, 'Cloudy'])/i}
    return(data.frame("Rain_relative_frequency" = relative_freq_Rain, "Cloudy_relative_frequency" = relative_freq_Cloudy))}


#2&3. Implement the Gibbs sampler sketched above for the Bayesian network in Figure 1 and draw
#100 samples from the joint probability distribution P(R, C | S = T, W = T).
#Estimate the marginal probability of rain, given that the sprinkler is on and the grass is wet
#P(R = T | S = T, W = T) from the 100 samples. (1 point)
    drawn_100_samples <- gibbs_sampler(100)
    cat("The marginal probability of P(R = T | S = T, W = T) calculated by drawing 100 samples from the joint 
            probability distribution P(R, C | S = T, W = T) = ", sum(drawn_100_samples['Rain']/nrow(drawn_100_samples)),"\n")
    
    
#4. Now draw 50,000 samples instead of 100 using the Gibbs sampler.
    sample_50000_1 = gibbs_sampler(50000)
    sample_50000_2 = gibbs_sampler(50000)
    cat("The marginal probability of P(R = T | S = T, W = T) calculated by drawing 50 000 samples from the joint 
            probability distribution P(R, C | S = T, W = T) = ", sum(sample_50000_1['Rain']/ nrow(sample_50000_1)),"\n")
    relative_freq_df_1 = count_relative_freq(sample_50000_1)
    relative_freq_df_2 = count_relative_freq(sample_50000_2)

    
#5. Provide the plot of the relative frequencies of R = T and C = T up to each iteration t against
#t, for two independent runs of the sampler. Suggest a burn-in time based on this plot. (1point)
    plot_of_rain_1 = ggplot(relative_freq_df_1,aes(x=(1:nrow(relative_freq_df_1)), y=Rain_relative_frequency)) +
        geom_line(stat = "identity", color = "blue") +
        ylim(0,1) +
        labs(title = "Run no.1 of sampler (with respect to Rain)",
             y = "Realative frequency of R = T",
             x = "Number of iteration")
    
    plot_of_rain_2 = ggplot(relative_freq_df_2,aes(x=(1:nrow(relative_freq_df_2)), y=Rain_relative_frequency)) +
        geom_line(stat = "identity", color = "blue") +
        ylim(0,1) +
        labs(title = "Run no.2 of sampler (with respect to Rain)",
             y = "Realative frequency of R = T",
             x = "Number of iteration")
    
    
    plot_of_cloudy_1 = ggplot(relative_freq_df_1,aes(x=(1:nrow(relative_freq_df_1)), y=Cloudy_relative_frequency)) +
        geom_line(stat = "identity", color = "red") +
        ylim(0,1) +
        labs(title = "Run no.1 of sampler (with respect to Cloudy)",
             y = "Realative frequency of C = T",
             x = "Number of iteration")
    
    plot_of_cloudy_2 = ggplot(relative_freq_df_2,aes(x=(1:nrow(relative_freq_df_2)), y=Cloudy_relative_frequency)) +
        geom_line(stat = "identity", color = "red") +
        ylim(0,1) +
        labs(title = "Run no.2 of sampler (with respect to Cloudy)",
             y = "Realative frequency of C = T",
             x = "Number of iteration")
    
    grid.arrange(plot_of_rain_1, plot_of_rain_2, plot_of_cloudy_1, plot_of_cloudy_2, nrow=2)

    
#6. Apply the Gelman test and plot potential scale reduction factor changes over the iterations
#using gelman.plot() from the coda package. Roughly speaking, this factor measures the ratio
#of the variances within and between independent runs of the sampler. Thus, for a stationary
#distribution, this factor should be close to 1.0. Suggest a burn-in time based on this plot. (2 points)
    mcmc_list = mcmc.list(mcmc(sample_50000_1), mcmc(sample_50000_2))
    gelman.plot(mcmc_list, autoburnin = FALSE)

    
#7. Investigate the auto-correlation among the samples. We expect adjacent members from a
#Gibbs sampling sequence to be positively correlated, and we can quantify the amount of
#this correlation by using the auto-correlation function. The lag-k auto-correlation ρk is the
#correlation between every draw and its kth neighbouring samples. Use the R-function acf() to
#provide plots for both variables Rain and Cloudy. Suggest an interval for drawing approximately
#independent samples. (2 points)
    acf(mcmc(sample_50000_1))
    acf(mcmc(sample_50000_2))

    
#8. Re-estimate P(R = T | S = T, W = T) based on 100 samples obtained after the suggested
#burn-in time and thinning-out. Compare with (3) and comment on your results. (1 point) 
    burn_in = 7000
    thinning_out = 4
    
    re_estimate = function(burn_in, thin_out){
        sample_size = burn_in + 100*thin_out
        samples_all = gibbs_sampler(sample_size)
        return(samples_all[seq(burn_in + thin_out, nrow(samples_all), thin_out), ])}
    
    re_estimated_results = re_estimate(as.integer(burn_in), as.integer(thinning_out))
    
    cat("Re-estimated P(R = T | S = T, W = T) based on 100 samples obtained after the suggested
    burn-in time and thinning-out. = ", sum(re_estimated_results['Rain']/nrow(re_estimated_results)),"\n")
