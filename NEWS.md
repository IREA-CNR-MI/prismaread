# prismaread 0.2.0 - 2019-06-04

* Major changes

    - prismaread now applies a GLT-based georeferencing to L1, L2B and L2C data, 
      based on the technique described in https://www.harrisgeospatial.com/docs/backgroundgltbowtiecorrection.html
      
    - Added functionality to create also outputs about acquisition angles and Lat/Lon;
    
    - Added the `prisma_extract_spectra` function to facilitate extraction of PRISMA
      data over features of vector files;
      
    - Added possibility to extract only some bands, using the `selbands_vnir` and 
      `selbands_swir` arguments to `convert_prisma`;
      
    - Removed the `out_file` argument, in favor of the `out_folder` and `out_basefile`
      arguments (see documentation); 
      
* Other changes

    - Miscellaneous bug fixes and optimization
    - Added pkgdown website

# prismaread 0.1.0 - 2019-03-25

* First beta "stable" version with several improvements

    - Proper conversion of L2D/L2C data, with results as reflectances after applying
      rescaling from the 0-65535 range
      
    - Proper conversion of L2B data, with results as radiances after applying
      rescaling from the 0-65535 range
      
    - Automatic application of scale and offset on L1 data, with results as radiances
    
    - Added argument `base_georef` to turn "ballpark" georeferencing of
      L1, L2B and L2C data on or off. Also, improved results thanks to proper chacking
      of corners
      
    - Check for proper "orientation" of imported L1, L2B and L2C data
    
    - Avoid saturation problems on L1 data by changing output data type

* Added a `NEWS.md` file to track changes to the package.

# prismaread 0.0.2 - 2019-10-02

*  First draft version


