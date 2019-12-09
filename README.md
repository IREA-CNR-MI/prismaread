
<!-- README.md is generated from README.Rmd. Please edit that file -->

# prismaread

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental)
[![Travis build
status](https://travis-ci.org/lbusett/prismaread.svg?branch=master)](https://travis-ci.org/lbusett/prismaread)
<!-- badges: end -->

The goal of prismaread is allowing to easily import PRISMA L1
hyperspectral images and convert them in a easier to use format (ENVI or
GeoTiff).

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("lbusett/prismaread")
```

## Work on Level 1 Data

The function to be used is `convert_prisma`. It takes as input the full
path of a PRISMA hdf5 image, an output file name and format, and a
series of switches allowing to decide which datasets should be saved. In
particular, the `FULL` argument allow to decide if a complete VNIR+SWIR
cube has to be created alongside the “single” VNIR and SWIR ones.

A “ballpark” georeferencing is provided, but it is not to be
trusted/used.

For example the following code:

``` r
library(prismaread)

in_file  <- "/home/lb/tmp/test/PRS_L1_STD_OFFL_20190825103112_20190825103117_0001.he5"
out_file <- "/home/lb/tmp/test/test_1"
out_format <- "ENVI"

# Save a full image, prioritizing the VNIR spectrometer and save in EVI format
convert_prisma(in_file    = in_file,
               out_file   = out_file,
               out_format = out_format,
               join_priority = "VNIR", 
               FULL = TRUE,
               PAN  = TRUE,
               CLOUD  = TRUE)
```

accesses the input file and saves both the VNIR and SWIR cubes, as well
as a full hyperspectral cube and the PAN and CLOUD images. **See
[documentation](reference/convert_prisma.html) of the convert\_prisma()
function for info on available arguments**.

The function also saves ancillary data related to wavelengths and fwhms
of the different images, and to hour and sun geometry at acquisition in
ancillary txt files.

### Creation of ATCOR files

Starting `v0.02` the function also allows automatic creation of text
files required to run an ATCOR atmospheric correction. Those files are
saved in the `ATCOR` subfolder of the main output folder.

In “standard” behaviour, only the three “standard” ATCOR files (`.wvl`,
`.dat` and `.cal`) are created, with the `.wvl` file containing nominal
wavelengths and FWHMs derived from the `cw` and `fwhm` attributes of the
*.he5* file- The user can however also choose to generate additional
ATCOR files, containing data about wavelengths and FWHMs related to
different “columns” of the data cube, as derived from the
`KDP_AUX/Cw_Vnir_Matrix`, `KDP_AUX/Cw_Swir_Matrix`,
`KDP_AUX/Cw_Fwhm_Matrix`, `KDP_AUX/Cw_Fwhm_Matrix` HDF layers. This
could allow running different atmospheric corrections for different
columns of the data, potentially allowing compensating “smile” effects
on the retrieved surface reflectances. For example:

``` r
library(prismaread)

in_file  <- "/home/lb/tmp/test/PRS_L1_STD_OFFL_20190825103112_20190825103117_0001.he5"
out_file <- "/home/lb/tmp/test/test_1"
out_format <- "ENVI"

# Save a full image, prioritizing the VNIR spectrometer and save in EVI format
convert_prisma(in_file    = in_file,
               out_file   = out_file,
               out_format = out_format,
               join_priority = "VNIR", 
               ATCOR = TRUE, 
               ATCOR_wls = c(200,800), 
               FULL = TRUE,
               PAN  = TRUE,
               CLOUD  = TRUE)
```

## Work on Level 2 Data

The function to be used is still `convert_prisma`, with similar syntax.
Only differences are:

1.  The “PAN”, “CLOUD”, “LC” and “GLINT” arguments are ignored, because
    they are not available in L2D files;
2.  The “ATCOR” and “ATCOR\_wls” arguments are ignored, because they are
    useless in this case;
3.  The output imagery is properly georeferenced.
