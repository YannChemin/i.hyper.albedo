## DESCRIPTION

*i.hyper.albedo* calculates broadband albedo from hyperspectral imagery
imported as 3D raster maps (`raster_3d`) by
[i.hyper.import](i.hyper.import.html).

The module reads wavelength metadata from hyperspectral 3D raster bands
and computes a weighted spectral integration to produce a single
broadband albedo raster map. Albedo represents the fraction of solar
radiation reflected by a surface, integrated across the solar spectrum
from 0.3 to 3.0 micrometers (300-3000 nm).

*i.hyper.albedo* is part of the **i.hyper** module family designed for
hyperspectral data import, processing, and analysis in GRASS. It
provides physically-based surface reflectance products for climate,
energy balance, and land surface studies.

The module implements spectral distance-weighted integration that
accounts for the spacing between successive spectral bands. This ensures
accurate representation of the continuous solar spectrum even when
sampling is irregular or sparse in certain wavelength regions.

Three weighting methods are available for spectral integration:

- **trapezoidal** -- Standard trapezoidal rule integration based on
  spectral bandwidth between adjacent bands (default)
- **solar** -- Solar irradiance-weighted integration using AM1.5 global
  spectrum approximation or custom solar spectrum
- **uniform** -- Equal weighting for all bands (simple arithmetic mean)

The trapezoidal method is recommended for general use as it properly
accounts for variable spectral sampling. Solar weighting provides the
most physically accurate albedo for energy balance applications, as it
weights each band by the amount of solar energy in that spectral region.

The output is a single-band raster map with values typically ranging
from 0 (complete absorption) to 1 (perfect reflection), though values
may exceed 1 in cases of high reflectance or sensor calibration
artifacts.

## NOTES

The module expects input data to be a 3D raster map created by
*i.hyper.import* or any 3D raster with wavelength metadata stored in
band-level metadata following the *i.hyper* standard format:
**wavelength**, **FWHM**, **valid**, and **unit**.

When wavelength metadata is not found, the module will fail with an
error message. It is essential to use properly imported hyperspectral
data with wavelength information.

The default wavelength range (300-3000 nm) covers the solar spectrum
relevant for albedo calculation: ultraviolet (300-400 nm), visible
(400-700 nm), near-infrared (700-1400 nm), and shortwave infrared
(1400-3000 nm). This range can be adjusted using the *min_wavelength*
and *max_wavelength* parameters, though values outside this range will
generate a warning.

The **-n** flag restricts processing to only bands marked as valid
(`valid=1`) in the metadata. This is useful for excluding bands known to
have quality issues, atmospheric absorption features, or other problems.

The **-i** flag provides information mode that prints band information
and computed weights without creating the output raster. This is helpful
for previewing which bands will be used and verifying the weighting
scheme before processing.

For solar-weighted albedo, the module uses a simplified AM1.5 global
solar spectrum approximation by default. For more accurate results,
users can provide a custom solar spectrum file in CSV format with two
columns: wavelength (nm) and irradiance (W/m²/nm). Standard solar
spectra are available from sources such as ASTM G173-03 or NREL.

The spectral distance calculation uses the trapezoidal rule: for each
band, the effective bandwidth is computed as half the distance to each
neighboring band. Edge bands use only the distance to their single
neighbor. This approach ensures proper integration even with irregular
spectral sampling.

Albedo values are stored as floating-point reflectance units (0-1
range). For visualization or further processing, users may want to scale
or classify the output using standard GRASS raster tools.

## EXAMPLES

::: code

    # Calculate broadband albedo from PRISMA data using trapezoidal integration
    i.hyper.albedo input=prisma \
                   output=prisma_albedo

    # Console output:
    Processing hyperspectral albedo for: prisma
    Wavelength range: 300 - 3000 nm (0.30 - 3.00 µm)
    Found 234 wavelength bands
    Bands in range: 234
    ======================================================================
    Band Information and Weights:
    ======================================================================
      Band | Wavelength (nm) |     FWHM |     Weight
    ----------------------------------------------------------------------
         1 |           402.71 |      8.5 | 0.003891
         2 |           412.65 |      8.4 | 0.003845
       ...
       234 |          2498.23 |     11.2 | 0.004123
    ======================================================================
    Calculating albedo using weighted integration...
    Successfully created albedo map: prisma_albedo
    --------------------------------------------------
    Albedo Statistics:
      Mean:   0.2847
      Min:    0.0312
      Max:    0.8923
      StdDev: 0.1245
    --------------------------------------------------
:::

::: code

    # Calculate solar-weighted albedo for energy balance studies
    i.hyper.albedo input=enmap \
                   output=enmap_albedo_solar \
                   weighting=solar

    # This weights bands by solar irradiance, giving more importance
    # to wavelengths where the sun emits more energy
:::

::: code

    # Calculate albedo only from visible spectrum (400-700 nm)
    i.hyper.albedo input=tanager \
                   output=tanager_visible_albedo \
                   min_wavelength=400 \
                   max_wavelength=700 \
                   weighting=trapezoidal
:::

::: code

    # Preview bands and weights without processing
    i.hyper.albedo input=prisma \
                   output=test \
                   weighting=solar \
                   -i

    # Info mode: displays band information and weights but creates no output
:::

::: code

    # Use only valid bands and custom wavelength range for NIR albedo
    i.hyper.albedo input=enmap \
                   output=enmap_nir_albedo \
                   min_wavelength=700 \
                   max_wavelength=1400 \
                   weighting=trapezoidal \
                   -n
:::

::: code

    # Solar-weighted albedo with custom solar spectrum
    i.hyper.albedo input=prisma \
                   output=prisma_albedo_custom \
                   weighting=solar \
                   solar_spectrum=custom \
                   solar_file=/data/ASTM_G173_AM15.csv
:::

::: code

    # Complete workflow: import and calculate albedo
    # Step 1: Import hyperspectral data
    i.hyper.import input=/data/PRISMA.he5 \
                   product=prisma \
                   output=prisma

    # Step 2: Calculate broadband albedo
    i.hyper.albedo input=prisma \
                   output=prisma_albedo \
                   weighting=solar

    # Step 3: Visualize albedo
    r.colors map=prisma_albedo color=grey
    d.rast map=prisma_albedo

    # Step 4: Calculate statistics by land cover class
    r.univar map=prisma_albedo zones=landcover
:::

## SEE ALSO

[i.hyper.import](i.hyper.import.html),
[i.hyper.rgb](i.hyper.rgb.html),
[i.hyper.composite](i.hyper.composite.html),
[i.hyper.preproc](i.hyper.preproc.html),
[i.hyper.explore](i.hyper.explore.html),
[i.hyper.export](i.hyper.export.html),
[r.mapcalc](https://grass.osgeo.org/grass-stable/manuals/r.mapcalc.html),
[r.univar](https://grass.osgeo.org/grass-stable/manuals/r.univar.html),
[r3.support](https://grass.osgeo.org/grass-stable/manuals/r3.support.html),
[r.colors](https://grass.osgeo.org/grass-stable/manuals/r.colors.html)

## REFERENCES

- Liang, S. (2001). Narrowband to broadband conversions of land surface
  albedo I: Algorithms. *Remote Sensing of Environment*, 76(2),
  213-238.
- Schaaf, C. B., et al. (2002). First operational BRDF, albedo nadir
  reflectance products from MODIS. *Remote Sensing of Environment*,
  83(1-2), 135-148.
- ASTM G173-03 (2020). Standard Tables for Reference Solar Spectral
  Irradiances: Direct Normal and Hemispherical on 37° Tilted Surface.
  *ASTM International*.
- Qu, Y., et al. (2014). Direct-estimation algorithm for mapping daily
  land-surface broadband albedo from MODIS data. *IEEE Transactions on
  Geoscience and Remote Sensing*, 52(2), 907-919.

## AUTHORS

Created for the i.hyper module family

Based on work by Alen Mangafić and Tomaž Žagar, Geodetic Institute of
Slovenia
