### POSSA

Power simulation for sequential analysis and multiple hypotheses.

This package is already in working condition as it is, but some parts of the underlying mechanism (and the related parameters) are still being optimized, hence future updates might bring some differences. The first definitive version will be published via CRAN (hopefully soon).

### Installation in R

Currently, POSSA is only available via GitHub. It can be installed as follows.

```
install.packages("devtools") # if "devtools" package is not yet installed
library("devtools")
install_github("gasparl/possa")
```

### Usage

Please see the vignettes:

- [General introduction](https://gasparl.github.io/possa/vignettes/intro.html "POSSA: Introduction")
- [Correlated samples](https://gasparl.github.io/possa/vignettes/correlations.html "POSSA: Correlated samples")
- [Multiple hypotheses](https://gasparl.github.io/possa/vignettes/multiple_hypotheses.html "POSSA: Multiple hypotheses")
- [Benchmarking](https://gasparl.github.io/possa/vignettes/benchmarking.html "POSSA: Benchmarking")

Upcoming:
- Real use examples (unequal samples and SDs; ranked data; DeLong's test; AOV)


For details about each function, see [the manual](https://github.com/gasparl/possa/blob/master/POSSA.pdf "POSSA manual") (or enter `help(xy)` or `?xy` in R for any specific function).


### Support

* If you run into an error despite carefully following the [documentation](https://github.com/gasparl/possa/blob/master/POSSA.pdf "POSSA.pdf"), [open a new issue](https://github.com/gasparl/possa/issues "Issues") or [email me](mailto:lkcsgaspar@gmail.com).
* If you have sound reason to believe that some of the presented statistics (or functions) are really not optimal and/or could be improved in some plausible way, [email me](mailto:lkcsgaspar@gmail.com).

### Citation

See `citation("POSSA")`.

[![R-CMD-check](https://github.com/gasparl/possa/workflows/R-CMD-check/badge.svg)](https://github.com/gasparl/possa/actions)
