### POSSA

Power simulation for sequential analysis and multiple hypotheses.

This package is already in working condition as it is, but some parts of the underlying mechanism (and the related parameters) are still being optimized, hence future updates might bring some differences. The first definitive version will be published via CRAN (hopefully soon).

### Installation in R

For now, POSSA is only available via GitHub. It can be installed as follows.

```
install.packages("devtools") # if "devtools" package is not yet installed
library("devtools")
install_github("gasparl/possa")
```

### Usage

Please see the vignettes:

1. [General introduction](https://gasparl.github.io/possa/vignettes/intro.html "POSSA: Introduction")
2. [Multiple hypotheses](https://gasparl.github.io/possa/vignettes/multiple_hypotheses.html "POSSA: Multiple hypotheses")
3. [Practical examples](https://gasparl.github.io/possa/vignettes/examples.html "POSSA: Practical examples") (unequal variances; ranked data; DeLong's test; AOV)
4. (Less important: [Benchmarking](https://gasparl.github.io/possa/vignettes/benchmarking.html "POSSA: Benchmarking"))

For detailed information about each function and parameter, see [the manual](https://github.com/gasparl/possa/blob/master/POSSA.pdf "POSSA manual").

### Support

* If you run into an error despite carefully following the [documentation](https://github.com/gasparl/possa/blob/master/POSSA.pdf "POSSA.pdf"), [open a new issue](https://github.com/gasparl/possa/issues "Issues").
* If you have sound reason to believe that some of the presented statistics or functions are not optimal and/or could be improved in some plausible way, [email me](mailto:lkcsgaspar@gmail.com).

### Contributions

Given that this project is not part of my "day job", the [issues](https://github.com/gasparl/possa/issues "Issues") grew beyond my capacity. Help would be most welcome; if you are up to it, feel free to [email me](mailto:lkcsgaspar@gmail.com) (or to reply to the relevant issue or open a new one).

### Citation

See `citation("POSSA")`.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6523855.svg)](https://doi.org/10.5281/zenodo.6523855) [![R-CMD-check](https://github.com/gasparl/possa/workflows/R-CMD-check/badge.svg)](https://github.com/gasparl/possa/actions)
