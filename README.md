# covid19-austria-omicron

Collection of various models for estimating the growth advantage and 
spread of the Omicron (B.1.1.529) COVID-19 variant in Austria. 

## Bayesian Multinomial GLMM

A Bayesian multinomial GLMM is fitted to the variant of concern data from Austria 
for estimating the growth advantage of Omicron (B.1.1.529) over Delta (B.1.617.2).
Currently the Delta, Omicron and Alpha variants of concern are included in the GLMM, 
since cases of those variants were reported in the investigated period. 

A weekly random effect is used in the GLMM to account for overdispersion in 
weekly case data. 

For reasoning behind this choice please see:
[Harrison, X. A. Using observation-level random effects to model overdispersion in count data in ecology and evolution. PeerJ 2, e616 (2014).](https://peerj.com/articles/616/).

The Bayesian multinomial GLMM is implemented using [brms](https://github.com/paul-buerkner/brms).

Model sources can be found in [/R/functions/model/multinomial_models.R](/R/functions/model/multinomial_models.R).

Fit of the multionmial GLMM to variants in Austria. 
![Multinomial GLMM fit](/output/austria-variants-multinomial-zoom.png)

Projection based upon the multinomial GLMM.
![Multinomial GLMM projection](/output/austria-variants-multinomial-projection.png)

Since the representativeness of sampling of variant of concern cases is not 
guaranteed we also estimated a probability of Omicron becoming dominant (more than half of all cases)
based upon a sampling adjustment. This adjustment is 0 for representative sampling,
meaning that none of the unsampled/unassigned cases are assumed to actually be Delta. 
For a sampling adjustment of 1 fully targeted sampling of Omicron cases is assumed, where 
all Omicron cases would have been found and all unsampled/unassigned cases are 
assumed to be Delta variant cases. 

Based upon this the posterior probability of Omicron being more than half of all cases is 
calculated. 
![Multinomial GLMM projection](/output/austria-more-than-half-omicron-sample_adj.png)

## Epidemia Time-Varying Effective Reproduction-Number Model

An experimental [epidemia](https://github.com/ImperialCollegeLondon/epidemia/) model
is also included, for estimating a time-varying reproduction number based upon a 
sampling factor introduced for accounting for the ratio of cases assigned to any 
variant of concern to all cases in each week. 

Based upon this samples are drawn from the posterior of rt for estimating a 
multiplicative advantage between the effective time varying reproduction numbers 
of the Delta and Omicron variants. 

This model is still very experimental and under development. 

## Quasibinomial GLM 

A [quasibinomial](https://rdrr.io/r/stats/family.html) GLM was also implemented 
to provide an additional cross-check for the Bayesian multinomial GLMM results. 

Fit of the quasibinomial GLM 
![Binomial GLM fit](/output/austria-variants-binomial-glm.png)

Fit shown on a log-odds scale
![Binomial GLM fit, log odds](/output/austria-variants-binomial-glm-log-odds.png)

## Variant plots

Various variant plots were also implemented.

An area plot showing the variant cases as a proportion of all cases. 
![Variants assigned](/output/austria-variants-assigned.png)

An area plot showing th variant cases as a proportion only of the cases assigned
to any variant of concern. 
![Variants only](/output/austria-variants-area.png)

## Build 

This project uses the [renv](https://rstudio.github.io/renv/articles/renv.html) 
and [targets](https://books.ropensci.org/targets/) R packages for reproducible
research. 

First you need to install [renv](https://rstudio.github.io/renv/articles/renv.html) 
and then install all required R packages using
```r
renv::restore()
```

The whole project can then be built using [targets](https://books.ropensci.org/targets/)
by just calling 
```r
tar_make()
```

Full locally parallelized build using targets futures is also supported using 
for example
```r
tar_make_future(workers=8)
```

## Data Sources

Variant of concern cases identified using variant specific PCR or sequencing 
are loaded from the weekly case data published by AGES.

[SARS-CoV-2-Varianten in Österreich](https://www.ages.at/themen/krankheitserreger/coronavirus/sars-cov-2-varianten-in-oesterreich/)
