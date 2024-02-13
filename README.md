# cSEAdb
Super Enhancer Fingerprints at the Constituent Level Across Cancers


## Installation

cSEAdb is an R package with its source code documented in this [repository](https://github.com/tenglab/cSEAdb).


```R
# Install cSEAdb
library(devtools)
devtools::install_github("https://github.com/tenglab/cSEAdb.git")
```

Or, Install with vignettes and dependencies.

```R
# Install cSEAdb with vignettes and dependencies
devtools::install_github("https://github.com/tenglab/cSEAdb.git",build_vignettes = TRUE)
```

## Using cSEAdb
First load cSEAdb,
```R
library(cSEAdb)
```
Then, follow the [**User's Guide**](https://github.com/tenglab/cSEAdb/blob/master/vignettes/cSEAdb.html) 
to load cSEAdb. In the guide, we detail how to query through cSEAdb
and to visualize SE signatures across cancers.

The users are also encouraged to refer to the help pages of R functions in this package. 
## Citation
If you use cSEAdb, please cite our paper at (https://doi.org/10.1371/journal.pcbi.1011873).

"Liu, Xiang, Nancy Gillis, Chang Jiang, Anthony McCofie, Timothy I. Shaw, Aik-Choon Tan, Bo Zhao, Lixin Wan, Derek R. Duckett, and Mingxiang Teng. "An Epigenomic fingerprint of human cancers by landscape interrogation of super enhancers at the constituent level." PLoS computational biology 20, no. 2: e1011873. DOI: https://doi.org/10.1371/journal.pcbi.1011873"
## Help
Feel free to leave any questions and bugs at [GitHub issues](https://github.com/tenglab/cSEAdb/issues).
