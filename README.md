<!-- README.md is generated from README.Rmd. Please edit that file -->
R package for the Straightforward Filtering INdeX (SFINX)
---------------------------------------------------------

**Affinity purification-mass spectrometry is one of the most common techniques for the analysis of protein-protein interactions, but inferring bona fide interactions from the resulting data sets remains notoriously difficult. We introduce SFINX, a Straightforward Filtering INdeX that identifies true-positive protein interactions in a fast, user-friendly, and highly accurate way. SFINX outperforms alternative techniques on benchmark data sets and is also available via the Web interface at <http://sfinx.ugent.be/>.**

Context
-------

The analysis of protein-protein interactions enables scientists to connect genotypes with phenotypes and to answer fundamental biological questions or generate new hypotheses on the functions of proteins. In this field, affinity purification-mass spectrometry is a classical approach wherein a protein of interest (bait) containing an epitope tag is purified under conditions that preserve the protein complex to allow the identification of co-purifying proteins by mass spectrometry.

Several software approaches already exist to separate the false-positives from the true-positives in these protein-protein interaction data sets, but none of these approaches combines high accuracy, speed and user-friendliness without the need for the input of external data. Therefore, we developed the Straightforward Filtering INdeX (SFINX), which excels at all these points.

Access
------

Users can easily access SFINX via the Web site interface at <http://sfinx.ugent.be/> or via this package. This package also allows users to more easily integrate SFINX in their own R pipelines or on their own servers.

You can install the released version of the package from CRAN using:

``` r
install.packages("sfinx")
```

To use the sfinx package that you installed in your library, you also have to load it as follows:

``` r
library(sfinx)
```

You can perform the standard SFINX analysis by using the sfinx() function of this package.

``` r
sfinx(DataInputExampleFile, BaitIdentityExampleFile)
```

Who
---

Kevin Titeca

More information
----------------

Further information and examples can be found in the sfinx-vignette file and by accessing the document information as follows:

``` r
?sfinx

help(sfinx)
```

The SFINX algorithm and its interface were published in the Journal of Proteome Research on January 4, 2016.

SFINX: Straightforward Filtering Index for Affinity Purification-Mass Spectrometry Data Analysis. Kevin Titeca, Pieter Meysman, Kris Gevaert, Jan Tavernier, Kris Laukens, Lennart Martens, and Sven Eyckerman. Journal of Proteome Research 2016 15 (1), 332-338. DOI: 10.1021/acs.jproteome.5b00666.

If you have suggestions or questions that remain after reading the article, the manual and the object information, you can contact us at <sfinxinteractomics@gmail.com> .
