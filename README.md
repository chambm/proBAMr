# proBAMr

[![TravisCI Build Status](https://travis-ci.org/chambm/proBAMr.svg?branch=master)](https://travis-ci.org/chambm/proBAMr)
[![Coverage Status](https://codecov.io/gh/chambm/proBAMr/branch/master/graph/badge.svg)](https://codecov.io/gh/chambm/proBAMr)
[![Platforms version](http://bioconductor.org/shields/availability/3.4/proBAMr.svg)](http://bioconductor.org/packages/devel/bioc/html/proBAMr.html)

Package Homepage: http://bioconductor.org/packages/devel/bioc/html/proBAMr.html 
Bug Reports: https://support.bioconductor.org/p/new/post/?tag_val=proBAMr

The interpretation of proteomics data is significantly enhanced by genomic annotation. 
Using proBAMr, peptide-spectrum-matches can be easily reannotated using user-specified gene annotation 
schemes and assembled into both protein and gene identifications. Using the genome as a common reference, 
proBAMr facilitates seamless proteomics and proteogenomics data integration. ProBAM files can be readily 
visualized in genome browsers and thus bring proteomics data analysis to a general audience beyond the 
proteomics community.

## Installation

To install this package, start R and enter:

```R
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("proBAMr")
```

Alternatively, you can install the latest version from github using devtools:

```R
# install.packages("devtools")
devtools::install_github("chambm/proBAMr")
```

## Documentation

To view documentation for the version of this package installed in your system, start R and enter:
```R
browseVignettes("proBAMr")
```