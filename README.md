# normHDP
A Bayesian Hierarchical Dirichlet Model for clustering.

Single-cell RNA sequencing (scRNA-seq) is a powerful technology that allows researchers to understand gene expression patterns at the single-cell level. However, analyzing
scRNA-seq data is challenging due to issues and biases in data collection. In normHDP, we constructed an integrated Bayesian model that simultaneously addresses normalization,
imputation and batch effects, and also non-parametrically clusters cells into groups which could be either unique or shared across multiple datasets.

Specifically, the Hierarchical Dirichlet process (HDP) is used to discover clusters of cells with similar mean expression and dispersion patterns. Posterior inferences are made based on a Gibbs sampler.

A full demo of how to use the R package is included in [here](https://github.com/jinluliu550/normHDP/R/full.demo.R). Code for generating data used for simulation 1 to 3 are shown in [here](https://github.com/jinluliu550/normHDP/R/simulations.R). Running the MCMC algorithm with chain width and chain depth set to 100 takes approximately 7 days on 10 cores to run for data with approximately 5,000 genes and 10,000 cells across 2 datasets. Running the MCMC algorithm with fixed allocation takes approximately 5 days to run 10,000 iterations on 10 cores.

