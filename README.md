# Bayesian MCMC phylogenetics tutorial in R
Fabricia Nascimento, Mario dos Reis and Ziheng Yang

This repository containts two tutorials that guide the user through writing simple MCMC phylogenetics software to estimate the molecular distance (and the transition/trasversion ratio) for a pairwise sequence alignment under the Jukes and Cantor (1969) (and Kimura's 1980) substitution model. The tutorials introduce concepts such as burn-in, mixing, convergence, efficiency and autocorrelation of the MCMC chain.

Directory `JC69/` contains the MCMC tutorial to calculate the molecular distance under the JC69 model. The directory contains three files. File `mcmc.JCd.R` contains the main R code with exercises. File `mcmc.JCrt.R` contains the solution to exercise 7 in the previous file. File `BayesianMCMC-JC.pdf` contains a more detailed explanation of the theory used in the tutorial.

Directory `K80/` contains the MCMC tutorial to calculate the molecular distance and the ts/tv rate under the K80 model. The tutorial is similar to the JC69 one, but focusing on a two parameter MCMC instead. File `mcmc.K80.R` contains the main R code. File `BayesianMCMC-K80.pdf` contains a detailed step-by-step explanation of the R code. The pdf file is the same as the webpage at:

https://thednainus.wordpress.com/2017/03/03/tutorial-bayesian-mcmc-phylogenetics-using-r/

In the K80 tutorial the user will be able to reproduce the plots to appear in our forthcoming review on MCMC phylogenetics (details coming soon!).
