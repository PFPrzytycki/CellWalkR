<img src="examples/CellWalkR_Vignette_files/figure-markdown_github/cellwalker_icon.png" id="id" class="class" width="50" height="50" /> CellWalkR
================

Installation
------------

Install CellWalkR for R using devtools as follows:

``` r
$ R
> install.packages("devtools")
> devtools::install_github("PFPrzytycki/CellWalkR")
```

Usage
-----

For a guide to using CellWalkR, see the provided [vignette](examples/CellWalkR_Vignette.md).

If you use CellWalkR please cite our paper:
Przytycki, P.F., Pollard, K.S. CellWalker integrates single-cell and bulk data to resolve regulatory elements across cell types in complex tissues. Genome Biol 22, 61 (2021). <https://doi.org/10.1186/s13059-021-02279-1>

AWS + TensorFlow
----------------

CellWalkR can also be run on AWS which vastly simplifies the process of running on GPUs using TensorFlow. Using GPUs allows the code to run more than 15 times faster. For a guide to running CellWalkR on AWS using GPUs see this [vignette](examples/CellWalkR_TensorFlow_Vignette.md).
