### POSSA: Power simulation for sequential analyses and multiple hypotheses

This R package serves to calculate, via simulation, power and appropriate "stopping" alpha boundaries (and/or futility bounds) for sequential analyses (i.e., "group sequential design") as well as for multiple hypotheses (multiple tests included in an analysis), given any specified global error rate. This enables the sequential use of practically any significance test, as long as the underlying data can be simulated in advance to a reasonable approximation.

For how to use, see the vignettes linked below.

### Installation in R

To install the stable version from [CRAN](https://cran.r-project.org/package=POSSA "The Comprehensive R Archive Network"), just run:

```R
install.packages("POSSA")
```

Then load with: `library("POSSA")`

Alternatively, if you want to install the latest (and potentially unstable) version from this repository:

```R
install.packages("devtools") # if "devtools" package is not yet installed
library("devtools")
install_github("gasparl/possa")
```

### Usage

Please see the vignettes:

1. [General introduction](https://gasparl.github.io/possa/vignettes/v_1_intro.html "POSSA: Introduction")
2. [Multiple hypotheses](https://gasparl.github.io/possa/vignettes/v_2_multiple_hypotheses.html "POSSA: Multiple hypotheses")
3. [Practical examples](https://gasparl.github.io/possa/vignettes/v_3_examples.html "POSSA: Practical examples") (unequal variances; ranked data; ANOVA; DeLong's test)

(Less important: [Benchmarking](https://gasparl.github.io/possa/vignettes/v_4_benchmarking.html "POSSA: Benchmarking"))

For detailed information about each function and parameter, see [the manual](https://github.com/gasparl/possa/blob/master/POSSA.pdf "POSSA manual").

### Support

* If you run into an error despite carefully following the [documentation](https://github.com/gasparl/possa/blob/master/POSSA.pdf "POSSA.pdf"), [open a new issue](https://github.com/gasparl/possa/issues "Issues").
* If you have sound reason to believe that some of the presented statistics or functions are not optimal and/or could be improved in some plausible way, [email me](mailto:lkcsgaspar@gmail.com).

### Contributions

Given that this project is not part of my "day job", the [intended additional features](https://github.com/gasparl/possa/issues "Issues") grew beyond my capacity. Help would be most welcome; if you are up to it, feel free to [email me](mailto:lkcsgaspar@gmail.com) (or to reply to the relevant issue or open a new one; see details [here](https://github.com/gasparl/possa/blob/master/CONTRIBUTING.md "CONTRIBUTING")).

### Etymology

The name "POSSA" refers not to the Italian meaning ["power, strength" (/ˈpɔs.sa/)](https://en.wiktionary.org/wiki/possa#Italian), but to the plural of ["possum" (/ˈpɑsə/)](https://en.wiktionary.org/wiki/possa#English).

### Citation

When you use POSSA in a publication, you can either cite the specific version you used (enter `citation("POSSA")` in R), or the following paper:

Lukács, G. (2022). POSSA: Power simulation for sequential analyses and multiple hypotheses. _Journal of Open Source Software, 7_(76), 4643. https://doi.org/10.21105/joss.04643

[![DOI](https://joss.theoj.org/papers/10.21105/joss.04643/status.svg)](https://doi.org/10.21105/joss.04643) ![](https://www.r-pkg.org/badges/version-last-release/POSSA "POSSA CRAN last version") ![](http://cranlogs.r-pkg.org/badges/POSSA?color=8585ad "POSSA CRAN monthly download count") [![R-CMD-check](https://github.com/gasparl/possa/workflows/R-CMD-check/badge.svg)](https://github.com/gasparl/possa/actions) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6523855.svg)](https://doi.org/10.5281/zenodo.6523855) [![codecov](https://codecov.io/gh/gasparl/possa/branch/master/graph/badge.svg?token=YVA1OCIDD7)](https://app.codecov.io/gh/gasparl/possa)
