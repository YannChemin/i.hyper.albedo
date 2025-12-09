#!/usr/bin/env python
##############################################################################
# MODULE:    i.hyper.albedo
# AUTHOR(S): Created for hyperspectral albedo computation
# PURPOSE:   Calculate broadband albedo from hyperspectral imagery
# COPYRIGHT: (C) 2025 by the GRASS Development Team
# SPDX-License-Identifier: GPL-2.0-or-later
##############################################################################

# %module
# % description: Calculate broadband albedo from hyperspectral imagery using spectral distance-weighted integration
# % keyword: imagery
# % keyword: hyperspectral
# % keyword: albedo
# % keyword: reflectance
# %end

# %option G_OPT_R3_INPUT
# % key: input
# % required: yes
# % description: Input hyperspectral 3D raster map (from i.hyper.import)
# % guisection: Input
# %end

# %option G_OPT_R_OUTPUT
# % key: output
# % required: yes
# % description: Output albedo raster map
# % guisection: Output
# %end

# %option
# % key: min_wavelength
# % type: double
# % required: no
# % answer: 300
# % description: Minimum wavelength to include (nanometers)
# % guisection: Wavelength Range
# %end

# %option
# % key: max_wavelength
# % type: double
# % required: no
# % answer: 3000
# % description: Maximum wavelength to include (nanometers)
# % guisection: Wavelength Range
# %end

# %option
# % key: weighting
# % type: string
# % required: no
# % options: trapezoidal,solar,uniform
# % answer: trapezoidal
# % description: Weighting method for spectral integration
# % guisection: Processing
# %end

# %option
# % key: solar_spectrum
# % type: string
# % required: no
# % options: default,custom
# % answer: default
# % description: Solar irradiance spectrum to use (for solar weighting)
# % guisection: Processing
# %end

# %option G_OPT_F_INPUT
# % key: solar_file
# % required: no
# % description: Path to custom solar spectrum file (wavelength,irradiance CSV)
# % guisection: Processing
# %end

# %flag
# % key: n
# % description: Only include bands marked as valid (valid=1)
# % guisection: Processing
# %end

# %flag
# % key: i
# % description: Print information about bands and wavelengths without processing
# % guisection: Processing
# %end

import sys
import os
import grass.script as gs
import numpy as np


def get_raster3d_info(raster3d):
    """Get information about 3D raster"""
    try:
        info = gs.raster3d_info(raster3d)
        return info
    except Exception as e:
        gs.fatal(f"Cannot get info for 3D raster {raster3d}: {e}")


def get_band_metadata(raster3d, band_num):
    """Extract metadata for a specific band"""
    band_name = f"{raster3d}#{band_num}"
    try:
        # Get band metadata using r3.support or similar
        metadata = gs.parse_command('r.support', map=band_name, flags='m')
        return metadata
    except:
        return {}


def parse_wavelength_from_metadata(raster3d, band_num):
    """Parse wavelength and validity from band metadata"""
    band_name = f"{raster3d}#{band_num}"
    wavelength = None
    fwhm = None
    valid = True
    unit = "nm"
    
    try:
        # Try to read metadata
        result = gs.read_command('r.support', map=band_name, flags='n')
        
        for line in result.split('\n'):
            line = line.strip()
            if line.startswith('wavelength='):
                wavelength = float(line.split('=')[1])
            elif line.startswith('FWHM='):
                fwhm = float(line.split('=')[1])
            elif line.startswith('valid='):
                valid = int(line.split('=')[1]) == 1
            elif line.startswith('unit='):
                unit = line.split('=')[1].strip()
    except:
        pass
    
    return wavelength, fwhm, valid, unit


def convert_wavelength_to_nm(wavelength, unit):
    """Convert wavelength to nanometers"""
    unit = unit.lower().strip()
    
    if unit in ['nm', 'nanometer', 'nanometers']:
        return wavelength
    elif unit in ['um', 'µm', 'micrometer', 'micrometers', 'micron', 'microns']:
        return wavelength * 1000.0
    elif unit in ['m', 'meter', 'meters']:
        return wavelength * 1e9
    else:
        gs.warning(f"Unknown wavelength unit '{unit}', assuming nanometers")
        return wavelength


def get_all_band_wavelengths(raster3d, only_valid=False):
    """Extract all band wavelengths and metadata from 3D raster"""
    info = get_raster3d_info(raster3d)
    depths = int(info['depths'])
    
    bands = []
    
    gs.verbose(f"Scanning {depths} bands for wavelength metadata...")
    
    for i in range(1, depths + 1):
        wavelength, fwhm, valid, unit = parse_wavelength_from_metadata(raster3d, i)
        
        if wavelength is not None:
            wavelength_nm = convert_wavelength_to_nm(wavelength, unit)
            
            if only_valid and not valid:
                gs.verbose(f"Band {i}: {wavelength_nm} nm - SKIPPED (invalid)")
                continue
            
            bands.append({
                'band_num': i,
                'wavelength': wavelength_nm,
                'fwhm': fwhm if fwhm else 0,
                'valid': valid,
                'unit': unit
            })
            
            gs.verbose(f"Band {i}: {wavelength_nm} nm (FWHM: {fwhm}, valid: {valid})")
    
    if not bands:
        gs.fatal("No wavelength metadata found in 3D raster bands. "
                "Please use data imported with i.hyper.import or add wavelength metadata.")
    
    # Sort bands by wavelength
    bands.sort(key=lambda x: x['wavelength'])
    
    return bands


def filter_bands_by_wavelength(bands, min_wl, max_wl):
    """Filter bands within specified wavelength range"""
    filtered = [b for b in bands if min_wl <= b['wavelength'] <= max_wl]
    
    gs.message(f"Wavelength range: {min_wl} - {max_wl} nm")
    gs.message(f"Total bands available: {len(bands)}")
    gs.message(f"Bands in range: {len(filtered)}")
    
    if not filtered:
        gs.fatal(f"No bands found in wavelength range {min_wl} - {max_wl} nm")
    
    return filtered


def calculate_spectral_distances(bands):
    """Calculate spectral distances between successive bands for weighting"""
    if len(bands) < 2:
        return [1.0]
    
    distances = []
    
    for i in range(len(bands)):
        if i == 0:
            # First band: distance to next band
            dist = (bands[i+1]['wavelength'] - bands[i]['wavelength']) / 2.0
        elif i == len(bands) - 1:
            # Last band: distance from previous band
            dist = (bands[i]['wavelength'] - bands[i-1]['wavelength']) / 2.0
        else:
            # Middle bands: average of distances to neighbors (trapezoidal rule)
            dist = (bands[i+1]['wavelength'] - bands[i-1]['wavelength']) / 2.0
        
        distances.append(dist)
    
    return distances


def get_solar_irradiance(wavelength, solar_spectrum='default'):
    """Get solar irradiance at given wavelength (W/m²/nm)"""
    # Simplified AM1.5 solar spectrum approximation
    # Peak around 500 nm, decreases in UV and IR
    
    if solar_spectrum == 'default':
        # Gaussian-like approximation of solar spectrum
        peak_wl = 500.0  # Peak around 500 nm
        if wavelength < 300:
            return 0.0
        elif wavelength < 400:
            # UV region - lower irradiance
            return 0.8 * np.exp(-((wavelength - peak_wl) / 300.0) ** 2)
        elif wavelength < 1000:
            # Visible to near-IR - peak region
            return 1.5 * np.exp(-((wavelength - peak_wl) / 400.0) ** 2)
        elif wavelength < 2500:
            # SWIR - decreasing
            return 0.6 * np.exp(-((wavelength - 1000) / 800.0) ** 2)
        else:
            # Far IR - very low
            return 0.1 * np.exp(-((wavelength - 2000) / 500.0) ** 2)
    
    return 1.0  # Fallback uniform weighting


def load_custom_solar_spectrum(filepath):
    """Load custom solar spectrum from CSV file"""
    try:
        import csv
        spectrum = {}
        with open(filepath, 'r') as f:
            reader = csv.reader(f)
            next(reader, None)  # Skip header if present
            for row in reader:
                if len(row) >= 2:
                    wl = float(row[0])
                    irrad = float(row[1])
                    spectrum[wl] = irrad
        return spectrum
    except Exception as e:
        gs.fatal(f"Failed to load solar spectrum file: {e}")


def calculate_weights(bands, weighting_method, solar_spectrum='default', solar_file=None):
    """Calculate weights for each band based on spectral distances and optional solar weighting"""
    distances = calculate_spectral_distances(bands)
    weights = []
    
    if weighting_method == 'uniform':
        # Uniform weighting
        total = sum(distances)
        weights = [d / total for d in distances]
        
    elif weighting_method == 'trapezoidal':
        # Trapezoidal rule - weight by spectral bandwidth
        total = sum(distances)
        weights = [d / total for d in distances]
        
    elif weighting_method == 'solar':
        # Solar irradiance weighted
        custom_spectrum = None
        if solar_spectrum == 'custom' and solar_file:
            custom_spectrum = load_custom_solar_spectrum(solar_file)
        
        solar_weights = []
        for i, band in enumerate(bands):
            wl = band['wavelength']
            if custom_spectrum and wl in custom_spectrum:
                solar_weight = custom_spectrum[wl]
            else:
                solar_weight = get_solar_irradiance(wl, solar_spectrum)
            
            # Combine solar weight with spectral distance
            solar_weights.append(solar_weight * distances[i])
        
        total = sum(solar_weights)
        if total > 0:
            weights = [w / total for w in solar_weights]
        else:
            gs.fatal("Solar weighting resulted in zero total weight")
    
    return weights


def print_band_info(bands, weights):
    """Print information about bands and their weights"""
    gs.message("=" * 70)
    gs.message("Band Information and Weights:")
    gs.message("=" * 70)
    gs.message(f"{'Band':>6} | {'Wavelength (nm)':>16} | {'FWHM':>8} | {'Weight':>10}")
    gs.message("-" * 70)
    
    for i, (band, weight) in enumerate(zip(bands, weights)):
        fwhm_str = f"{band['fwhm']:.1f}" if band['fwhm'] else "N/A"
        gs.message(f"{band['band_num']:>6} | {band['wavelength']:>16.2f} | {fwhm_str:>8} | {weight:>10.6f}")
    
    gs.message("=" * 70)
    gs.message(f"Total bands: {len(bands)}")
    gs.message(f"Wavelength range: {bands[0]['wavelength']:.1f} - {bands[-1]['wavelength']:.1f} nm")
    gs.message(f"Total weight: {sum(weights):.6f}")
    gs.message("=" * 70)


def calculate_albedo(raster3d, bands, weights, output):
    """Calculate weighted albedo from bands"""
    gs.message("Calculating albedo using weighted integration...")
    
    # Build mapcalc expression for weighted sum
    terms = []
    for band, weight in zip(bands, weights):
        band_name = f"{raster3d}#{band['band_num']}"
        terms.append(f"({weight} * {band_name})")
    
    expression = f"{output} = " + " + ".join(terms)
    
    gs.verbose(f"Mapcalc expression length: {len(expression)} characters")
    gs.verbose(f"Number of terms: {len(terms)}")
    
    # Execute r.mapcalc
    try:
        gs.run_command('r.mapcalc', expression=expression, overwrite=True)
        gs.message(f"Successfully created albedo map: {output}")
    except Exception as e:
        gs.fatal(f"Failed to calculate albedo: {e}")
    
    # Set metadata
    gs.run_command('r.support', map=output, 
                  title="Broadband Albedo",
                  units="reflectance",
                  description=f"Spectral albedo integrated from {bands[0]['wavelength']:.0f} to {bands[-1]['wavelength']:.0f} nm")
    
    # Calculate and report statistics
    try:
        stats = gs.parse_command('r.univar', map=output, flags='g')
        gs.message("-" * 50)
        gs.message("Albedo Statistics:")
        gs.message(f"  Mean:   {float(stats['mean']):.4f}")
        gs.message(f"  Min:    {float(stats['min']):.4f}")
        gs.message(f"  Max:    {float(stats['max']):.4f}")
        gs.message(f"  StdDev: {float(stats['stddev']):.4f}")
        gs.message("-" * 50)
    except:
        pass


def main(options, flags):
    """Main function"""
    input_raster = options['input']
    output_raster = options['output']
    min_wl = float(options['min_wavelength'])
    max_wl = float(options['max_wavelength'])
    weighting = options['weighting']
    solar_spectrum = options['solar_spectrum']
    solar_file = options.get('solar_file', None)
    only_valid = flags['n']
    info_only = flags['i']
    
    # Validate wavelength range (0.3 to 3.0 micrometers = 300 to 3000 nm)
    if min_wl < 300 or max_wl > 3000:
        gs.warning(f"Wavelength range extends beyond typical albedo range (300-3000 nm)")
    
    if min_wl >= max_wl:
        gs.fatal("Minimum wavelength must be less than maximum wavelength")
    
    gs.message(f"Processing hyperspectral albedo for: {input_raster}")
    gs.message(f"Wavelength range: {min_wl} - {max_wl} nm ({min_wl/1000:.2f} - {max_wl/1000:.2f} µm)")
    
    # Get all bands with wavelength metadata
    all_bands = get_all_band_wavelengths(input_raster, only_valid=only_valid)
    
    # Filter bands by wavelength range
    bands = filter_bands_by_wavelength(all_bands, min_wl, max_wl)
    
    # Calculate weights
    weights = calculate_weights(bands, weighting, solar_spectrum, solar_file)
    
    # Print information
    print_band_info(bands, weights)
    
    if info_only:
        gs.message("Info mode: No albedo map created.")
        return 0
    
    # Calculate albedo
    calculate_albedo(input_raster, bands, weights, output_raster)
    
    gs.message(f"Albedo calculation complete: {output_raster}")
    gs.message(f"Method: {weighting} weighting")
    
    return 0


if __name__ == "__main__":
    options, flags = gs.parser()
    sys.exit(main(options, flags))
