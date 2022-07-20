### POSSA: Power simulation for sequential analysis and multiple hypotheses

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

1. [General introduction](https://gasparl.github.io/possa/vignettes/1_intro.html "POSSA: Introduction")
2. [Multiple hypotheses](https://gasparl.github.io/possa/vignettes/2_multiple_hypotheses.html "POSSA: Multiple hypotheses")
3. [Practical examples](https://gasparl.github.io/possa/vignettes/3_examples.html "POSSA: Practical examples") (unequal variances; ranked data; ANOVA; DeLong's test)

(Less important: [Benchmarking](https://gasparl.github.io/possa/vignettes/4_benchmarking.html "POSSA: Benchmarking"))

For detailed information about each function and parameter, see [the manual](https://github.com/gasparl/possa/blob/master/POSSA.pdf "POSSA manual").

### Support

* If you run into an error despite carefully following the [documentation](https://github.com/gasparl/possa/blob/master/POSSA.pdf "POSSA.pdf"), [open a new issue](https://github.com/gasparl/possa/issues "Issues").
* If you have sound reason to believe that some of the presented statistics or functions are not optimal and/or could be improved in some plausible way, [email me](mailto:lkcsgaspar@gmail.com).

### Contributions

Given that this project is not part of my "day job", the [intended additional features](https://github.com/gasparl/possa/issues "Issues") grew beyond my capacity. Help would be most welcome; if you are up to it, feel free to [email me](mailto:lkcsgaspar@gmail.com) (or to reply to the relevant issue or open a new one).

### Citation

See `citation("POSSA")`.

![](https://www.r-pkg.org/badges/version-last-release/POSSA "POSSA CRAN last version") [![R-CMD-check](https://github.com/gasparl/possa/workflows/R-CMD-check/badge.svg)](https://github.com/gasparl/possa/actions) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6523855.svg)](https://doi.org/10.5281/zenodo.6523855)
