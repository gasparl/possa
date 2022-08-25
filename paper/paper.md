---
title: 'POSSA: Power simulation for sequential analyses and multiple hypotheses'
tags:
  - r
  - statistical power
  - simulation
  - sequential design
  - multiple hypotheses
  - multiple testing
authors:
 - name: Gáspár Lukács
   orcid: 0000-0001-9401-4830
   affiliation: 1
affiliations:
 - name: Tilburg University, Department of Methodology and Statistics, The Netherlands
   index: 1
date: 20 July 2022
bibliography: paper.bib
---

# Summary

Conventional experiment designs have a fixed sample size: observations are collected until the planned sample size is reached, and then the experiment is concluded. Alternatively however, sequential designs can also be used [e.g., @lakens:2014]. In that case, interim analyses can be performed during data collection, and whenever the analysis suggests sufficient evidence for the presence or absence of the given effect, the experiment can be concluded without further data collection. This can drastically reduce the required sample sizes, saving time, expenses, and effort. However, sequential analyses require adjustments for Type 1 error rate (i.e., the ratio of false significant findings) and affect statistical power [which is crucial to determine the required sample size, e.g., @cohen:1988; @perugini:2018]. Furthermore, regardless of sequential analysis, in case of multiple hypotheses (multiple tests included in an analysis), for both Type 1 error rate and power, correction is necessary [e.g., @khandis:2020].

The present `POSSA` R package serves to perform, via a simulation framework, the necessary adjustments and calculate power for sequential analyses as well as for multiple hypotheses for practically any significance test.

# Statement of need

While sequential designs can be extremely useful, their application is limited by the lack of software that performs the necessary statistical adjustments. Recently, there have been several new related software solutions [e.g., the R packages by @anderson:2022; @pahl:2022; @wassmer:2022; for more, see @weigl:2020], but all of these apply only to single and specific parametric tests. Solutions for simulation-based power analysis also exist [e.g., @green:2016; @hughes:2017; @lakens:2021], but none specifically for sequential designs or multiple hypotheses. The `POSSA` R package allows power and Type 1 error calculations and corrections for analyses in fixed as well as sequential designs including any, and any number of, statistical null hypothesis tests.

# Example

A simplistic example is presented below to demonstrate how `POSSA` works.

First, create a custom sampling function, simulating the data that would be obtained in an actual single run of the experiment in case of (a) existing true effect in the population, and (b) no effect. In this example, it is simply two normally distributed sets of numbers.

```r
customSample = function(sampleSize) {
  list(
    # baseline
    sample1 = rnorm(sampleSize, mean = 0, sd = 10),
    # no effect
    sample2_h0 = rnorm(sampleSize, mean = 0, sd = 10),
    # true effect
    sample2_h1 = rnorm(sampleSize, mean = 5, sd = 10)
  )
}
```

Second, create a corresponding custom testing function, including all tests that would be conducted on the obtained samples. Here, it is simply a single _t_-test.

```r
customTest = function(sample1, sample2_h0, sample2_h1) {
  c(
    # no effect
    p_h0 = t.test(sample1, sample2_h0, 'less', var.equal = TRUE)$p.value,
    # true effect
    p_h1 = t.test(sample1, sample2_h1, 'less', var.equal = TRUE)$p.value
  )
}
```

Third, pass the sampling and testing functions to the `POSSA::sim` function, with any desired observation numbers (at the interim stops and at the final stop). This will simulate and return a great number of _p_ values.

```r
dfPvalsSeq = sim(fun_obs = customSample,
                n_obs = c(27, 54, 81),
                fun_test = customTest)
```

Fourth and last, pass the obtained `dfPvalsSeq` `data.frame` to the `POSSA::pow` function. In this case, `alpha_locals = NA` specifies that all local alphas should be identical [Pocock's correction, @pocock:1977].

```r
pow(dfPvalsSeq, alpha_locals = NA)
```

The information returned by this function includes the values of the local alphas (to achieve, by default, a Type 1 error rate of `.05`), the actual Type 1 error rate, the statistical power, and the average expected sample sizes in case of a true effect and in case of no effect.

The same workflow can be used to accommodate any other analyses.

# References
