# normHDP
A Bayesian Hierchical Dirichlet Model for clustering.

Single-cell RNA sequencing (scRNA-seq) is a powerful technology that allows researchers to understand gene expression patterns at the single-cell level. However, analysing
scRNA-seq data is challenging due to issues and biases in data collection. In normHDP, we constructed an integrated Bayesian model that simultaneously addresses normalization,
imputation and batch effects, and also nonparametrically clusters cells into groups which could be either unique or shared across multiple datasets.

Specifically, the Hierachical Dirichlet process (HDP) is used to discover clusters of cells with similar mean expression and dispersion patterns. Posterior inferences are
made based on a Gibbs sampler.
