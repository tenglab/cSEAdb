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

## Help
Feel free to leave any questions and bugs at [GitHub issues](https://github.com/tenglab/cSEAdb/issues).
