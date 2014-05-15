Healy.et.al-2014-Longevity-ProcB
================================

These following R functions were used in Healy et al 2014 Proc.B.
They allow to run comparative analysis for multiple trees from a
Bayesian tree building distribution (instead of using the consensus tree).

The code is based on MCMCglmm package.

"Rather than basing our analyses on just a single phylogenetic
tree and assuming this tree was known without error, we instead
used a distribution of trees.
[...]
All analyses were carried out in R v. 3.0.2 [35]. Maximum
longevity and body mass (and BMR, see below) were log10 transformed
to correct inherent skewness before being mean centred
and expressed in units of standard deviation.
We fitted our models using Bayesian phylogenetic mixed
models from the MCMCglmm package [36], to account for nonindependence
in species traits introduced as a result of common
ancestry [18]. MCMCglmm uses a Markov chain Monte Carlo
(MCMC)estimation approach and accounts for non-independence
among closely related species by including the phylogenetic
relationships among species as a random variable.We determined
the number of iterations, thinning and the burn-in period for each
model run across all trees using diagnostics in the coda package
[37] and we checked for convergence between model chains
using the Gelman-Rubin statistic, the potential scale reduction
factor (PSR), with all models required to have a PSR below 1.1
[38]. Following the recommendations of Hadfield [36], we used
an uninformative inverse-Wishart distribution (with variance, V,
set to 0.5 and belief parameter, nu, set to 0.002) and a parameter
expanded prior, with a half-Cauchy distribution (described by
the parameters V ¼ 0.5, nu ¼ 1, the prior mean alpha.mu ¼ 0,
and alpha.V ¼ 102, which represents the prior standard deviation
with a scale of 10), for the random factor to improve mixing and
decrease autocorrelation among iterations."
