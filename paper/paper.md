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

Conventional experiment designs have a fixed sample size: observations are collected until the planned sample size is reached, and then the experiment is concluded. Alternatively however, sequential designs can also be used [e.g., @lakens:2014]. In that case, interim analyses can be performed during data collection, and whenever the analysis suggests sufficient evidence for the presence or absence of the given effect, the experiment can be concluded without further data collection. This can drastically reduce the required sample sizes, saving time, expenses, and effort. However, sequential analyses require adjustments for Type 1 error rate (i.e., the ratio of false significant findings) and affect statistical power [which is crucial to determine the required sample size, e.g., @perugini:2018; @cohen:1988]. Furthermore, regardless of sequential analysis, in case of multiple hypotheses (multiple tests included in an analysis), for both Type 1 error rate and power, correction is necessary [though frequently ignored; e.g., @khandis:2020]. 

The present `POSSA` R package serves to make the necessary adjustments and calculate power for sequential analyses as well as for multiple hypotheses for practically any significance test.

# Statement of need

While sequential designs can be extremely useful, their application is limited by the lack of software that performs the necessary statistical adjustments. Recently, there have been several new related software solutions [e.g., the R packages by @wassmer:2022; @anderson:2022; @rahl:2022; for more, see @weigl:2020], but all of these apply only to single and specific parametric tests. Solutions for simulation-based power analysis also exist [e.g., @hughes:2017; @green:2016; @lakens:2021], but none specifically for sequential designs or multiple hypotheses. The `POSSA` R package allows power and Type 1 error calculations and corrections for analyses including any, and any number of, statistical null hypothesis tests.

# References