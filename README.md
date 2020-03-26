
<!-- README.md is generated from README.Rmd. Please edit that file -->

# prismaread

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/lbusett/prismaread.svg?branch=master)](https://travis-ci.org/lbusett/prismaread)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

The goal of prismaread is allowing to easily import PRISMA hyperspectral
data (<http://www.prisma-i.it/index.php/it/>) and convert them to a
easier to use format (ENVI or GeoTiff).

# Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("lbusett/prismaread")
library(prismaread)
```

# Work on Level 1 Data

The function to be used is `convert_prisma`. It takes as input the full
path of a PRISMA hdf5 image, an output file name and format, and a
series of switches allowing to decide which datasets should be saved.

In particular, the `FULL` argument allows deciding if a complete
VNIR+SWIR cube has to be created alongside the “single” VNIR and SWIR
ones. In that case, the ‘join\_priority’ keyword is used to decide if
keeping bands from the “VNIR” or the “SWIR” data cube in the wavelength
were they overlap.

A “ballpark” georeferencing is provided, but it is **not** to be
trusted/used. (this may change in the near future)

For example the following
code:

``` r
in_file  = "/home/lb/tmp/test/PRS_L1_STD_OFFL_20190825103112_20190825103117_0001.he5"
out_file = "/home/lb/tmp/test/test_1"
out_format = "ENVI"

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
on the retrieved surface reflectances. For
example:

``` r
in_file  = "/home/lb/tmp/test/PRS_L1_STD_OFFL_20190825103112_20190825103117_0001.he5"
out_file = "/home/lb/tmp/test/test_1"
out_format = "ENVI"

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

**IMPORTANT NOTE**

The latter functionality may be only appliable to “HRC” L1 data cubes.
We are currently investigating this - proceed with caution\!

# Work on Level 2 Data

The function to be used is still `convert_prisma`, with similar syntax.
Only differences are:

1.  The “PAN”, “CLOUD”, “LC” and “GLINT” arguments are ignored, because
    they are not available in L2 files;
2.  The “ATCOR” and “ATCOR\_wls” arguments are ignored, because they are
    useless in this case;
3.  The output imagery is properly georeferenced, **if level is 2D**.
    Otherwise, the same “ballpark” georeferencing as for L1 data is
    applied.

# Output Formats

  - Outputs are provided as rasters in **ENVI** or **GEOTIFF** format
    according to user’s choice.

  - Filenames are built starting from the output file name provided by
    the user, by adding appropriate suffixes. For example, if the user
    specified `out_file = "D:/myoutfolder/myoutfil"`, and `source =
    "HCO"`, the output file for the VNIR cube will be
    `D:/myoutfolder/myoutfil_HCO_VNIR.envi` (or
    `D:/myoutfolder/myoutfil_HCO_VNIR.tif`).

  - Measure units of the output hyperspectral data are as follows:
    
    | LEVEL |  Variable   |    Measure Units    |
    | :---: | :---------: | :-----------------: |
    |  L1   |  Radiance   | W / m^2 \* sr \* um |
    |  L2B  |  Radiance   | W / m^2 \* sr \* um |
    |  L2C  | Reflectance |  unitless (ratio)   |
    |  L2D  | Reflectance |  unitless (ratio)   |
    

  - If output format is “ENVI”, the wavelengths of the different bands
    for the hyperspectral cubes are properly written in the appropriate
    header (.hdr) file.

  - Irrespective from output format, info about wavelengths and fwhms of
    the hyperspectral cubes are saved in appropriate txt files. For
    example, if the output file is
    `D:/myoutfolder/myoutfil_HCO_VNIR.envi`, info about the wavelengths
    is saved in `D:/myoutfolder/myoutfil_HCO_VNIR_meta.txt`

  - Info about acquisition date and angles is saved in a dedicated txt
    file. For example, if output file is
    `D:/myoutfolder/myoutfil_VNIR.envi`, info about the angles is saved
    in `D:/myoutfolder/myoutfil_HCO_VNIR_meta.txt`

# Future Work

  - Test possibility to more properly georeference L1/L2B/L2C data using
    the curvilinear grids functionality in package `stars`

  - Implement possibility to apply masks base on the ERR\_MATRIX cubes

  - Clean up code
