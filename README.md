# TFactSR

# Installation
--------------
```R
install.packages("devtools")

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("qvalue")

devtools::install_github("afukushima/TFactSR", build_vignettes = TRUE)
```

# Tutorial
--------------
```{R}
browseVignettes("TFactSR")
```

or see https://afukushima.github.io/TFactSR/
