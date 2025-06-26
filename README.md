This repository contains code adapted from the methods described in [manuscript](https://doi.org/10.1371/journal.pgen.1011022) to extend the original analyses. 

## Running the code and generating the results

To run the code and generate the results, follow these steps:

- Submit the `generate_data_parallel.sh` script to generate the data using paralel array jobs. 
- Submit `boom.sh, bslmm.sh, and one-at-a-time.sh` to run the analysis using BoomSpikeSlab, Bayesian Sparse Linear Mixed Model, and univariate regression, respectively. 
- Use `true_positive_rate_mse.R` to calculate/plot the true positive rate and mean squared error of the estimates from the different methods.

