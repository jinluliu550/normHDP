# normHDP
A Bayesian Hierchical Dirichlet Model for clustering.

Single-cell RNA sequencing (scRNA-seq) is a powerful technology that allows researchers to understand gene expression patterns at the single-cell level. However, analysing
scRNA-seq data is challenging due to issues and biases in data collection. In normHDP, we constructed an integrated Bayesian model that simultaneously addresses normalization,
imputation and batch effects, and also nonparametrically clusters cells into groups which could be either unique or shared across multiple datasets.

Specifically, the Hierachical Dirichlet process (HDP) is used to discover clusters of cells with similar mean expression and dispersion patterns. Posterior inferences are made based on a Gibbs sampler.

A full demo of how to use the R package is included in [Simulation 1 example - pdf](https://github.com/jinluliu550/normHDP/blob/master/man/github_case1.pdf) which is equivalent to Simulation 1 in the article; a simple example with 2 datasets was used in this case and the overall algorithm takes approximately 30 minutes to run on 4 cores (html version of this example is provided [here](https://github.com/jinluliu550/normHDP/blob/master/man/github_case1.html), some of the example outputs are included in the [data folder](https://github.com/jinluliu550/normHDP/blob/master/data)). An additional example is included [here](https://github.com/jinluliu550/normHDP/blob/master/man/github_case3.pdf) which demos a further example with pre-classified marker
genes (Simulation 3) which takes approximately two and half hours to run.
