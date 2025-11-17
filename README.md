# Some stan demo code

## Pre-requisites

- R > 4.4 (you can use Python but since I don't use cmdstanpy I may not be able to help you)
- Windows: Rtools [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/) or if you are on Windows 11, you can use 

```cmd
winget.exe install --id "RProject.Rtools"
```

- StanHeaders + rstan: via CRAN

```r
# install.packages('StanHeaders')
install.packages('rstan') 
```

- cmdstanr: this is a newer version of Stan interface on R wrapping around cmdstan. It supports Pathfinder [https://arxiv.org/abs/2108.03782](https://arxiv.org/abs/2108.03782) and Laplace Approximation but less inference tools and not on CRAN.

```r
# we recommend running this in a fresh R session or restarting your current session
install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
```

- brms: very powerful and hackable Stan model builder, cover 80%+ of normal use cases. Support mostly everything in lme4, mgcv, gamm4, quantreg and more. It generates readable Stan code so you can even learn and improve from it (my way).

```r
install.packages('brms')
```

## Other packages you may interest

- rstanarm: pre-compiled Stan models that does a lot of things brms does. It doesn't need compliation for every new model but it is less flexible and does not support mulit-threading.

```r
install.packages('rstanarm')
```

- loo: comparing models by estimated log pointwise density (elpd) using PSIS sample (Vehtari). This can be accessed from brms and rstanarm but if you build your custom model you may need it.

- bayesplot, tidybayes, ggdist: very powerful packages to plot your MCMC samples as well as diagnostic metrics. Can be accessed from brms and rstanarm but if you build your own custom model you may need it.

- projpred: variable selection using projection and elpd.


