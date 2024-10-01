#!/usr/bin/env python3
'''
Created on 10/02/2024
@author: M Nergizci, some functions adapted from M Lazecky.  
'''
import numpy as np
import py_gamma as pg
import subprocess
import shutil
import time
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from shapely.geometry import Point
from matplotlib.path import Path
from lics_unwrap import *
#import framecare as fc
from scipy import ndimage 
from geocube.api.core import make_geocube
import geopandas as gpd
import pandas as pd
#import daz_lib_licsar as dl
#import LiCSAR_misc as misc
import sys
import rasterio
from rasterio.merge import merge
import os
#import LiCSquery as lq
from scipy.constants import speed_of_light
from scipy.signal import convolve2d, medfilt
from scipy.ndimage import sobel
from scipy.signal import medfilt2d
from scipy.ndimage import gaussian_filter
import dask.array as da
from dask_image import ndfilters
from scipy.ndimage import generic_filter
import scipy.ndimage.filters as ndfilters
from scipy.interpolate import interp1d
pi=np.pi


def extract_burst_overlaps(frame):
    # Read GeoJSON data
    data_temp = gpd.read_file(frame + '.geojson')

    # Change CRS to EPSG:4326
    data_temp = data_temp.to_crs(epsg=4326)

    # Extract subswath information
    if frame.startswith('00'):
        data_temp['swath'] = data_temp.burstID.str[4]
    elif frame.startswith('0'):
        data_temp['swath'] = data_temp.burstID.str[5]
    else:
        data_temp['swath'] = data_temp.burstID.str[6]

    # Divide frame into subswaths
    data_temp = data_temp.sort_values(by=['burstID']).reset_index(drop=True)
    sw1 = data_temp[data_temp.swath == '1']
    sw2 = data_temp[data_temp.swath == '2']
    sw3 = data_temp[data_temp.swath == '3']

    # Divide burst overlaps into odd and even numbers
    a1 = sw1.iloc[::2]
    b1 = sw1.iloc[1::2]
    a2 = sw2.iloc[::2]
    b2 = sw2.iloc[1::2]
    a3 = sw3.iloc[::2]
    b3 = sw3.iloc[1::2]

    # Find burst overlaps
    overlap_gdf1 = gpd.overlay(a1, b1, how='intersection')
    overlap_gdf2 = gpd.overlay(a2, b2, how='intersection')
    overlap_gdf3 = gpd.overlay(a3, b3, how='intersection')

    # Merge swath overlaps
    gpd_overlaps = gpd.GeoDataFrame(pd.concat([overlap_gdf1, overlap_gdf2, overlap_gdf3], ignore_index=True))

    return gpd_overlaps, overlap_gdf1, overlap_gdf2, overlap_gdf3


##############################################
def open_geotiff(path, fill_value=0):
    '''
    This code help open geotiff with gdal and remove nan to zero!
    '''
    try:
        bovl = gdal.Open(path, gdal.GA_ReadOnly)
        if bovl is None:
            raise Exception("Failed to open the GeoTIFF file.")

        band = bovl.GetRasterBand(1)
        bovl_data = band.ReadAsArray()

        # Replace NaN values with the specified fill_value
        #print(bovl_data.dtype)
        bovl_data[np.isnan(bovl_data)] = fill_value

        return bovl_data
    except Exception as e:
        print("Error:", e)
        return None

def export_to_tiff(output_filename, data_array, reference_tif_path):
    """
    Export a NumPy array to a GeoTIFF file using a reference TIFF for geospatial properties.
    
    Parameters:
    - output_filename: String, the name of the output GeoTIFF file.
    - data_array: NumPy array containing the data to be exported.
    - reference_tif_path: String, the file path of the reference GeoTIFF.
        
    Returns:
    None

    # Example usage:
    # output_filename = 'exported_data.tif'
    # data_array = your_numpy_array_here  # NumPy array you want to export
    # reference_tif_path = 'path/to/reference.tif'
    # export_to_tiff_with_ref(output_filename, data_array, reference_tif_path)
    """
    # Open the reference TIFF to read its spatial properties
    ref_ds = gdal.Open(reference_tif_path)
    geotransform = ref_ds.GetGeoTransform()
    projection = ref_ds.GetProjection()

    driver = gdal.GetDriverByName("GTiff")
    
    # Get the dimensions of the data array
    row, col = data_array.shape
    
    # Create the output GeoTIFF with Deflate compression
    outdata = driver.Create(
        output_filename, col, row, 1, gdal.GDT_Float32,
        options=["COMPRESS=DEFLATE"]  # Add Deflate compression
    )
    
    # Set the geotransform and projection from the reference TIFF
    outdata.SetGeoTransform(geotransform)
    outdata.SetProjection(projection)
    
    # Write data to the raster band
    outdata.GetRasterBand(1).WriteArray(data_array)
    
    # Flush the cache to disk to write changes
    outdata.FlushCache()
    
    # Cleanup
    ref_ds = None
    outdata = None


def check_file_exists(file_path, file_description):
    """Check if a file exists and notify the user."""
    if os.path.exists(file_path):
        print(f"{file_description} found at {file_path}")
        return True
    else:
        print(f"{file_description} not found at {file_path}")
        return False

def gradient_nr(data, deramp=True):
    """Calculates gradient of continuous data (not tested for phase)

    Args:
        xar (np.ndarray): A NumPy array, e.g. ifg['unw']
        deramp (bool): If True, it will remove the overall ramp

    Returns:
        np.ndarray
        
        gradient=calculate_gradient(azof,deramp=False)
        plt.figure(figsize=(10,10))
        plt.imshow(gradient, cmap='viridis', vmax=0.5)
        plt.colorbar()    
    """
    gradis = data.copy()  # Assuming xar is already a NumPy array
    vgrad = np.gradient(gradis)  # Use NumPy's gradient function
    gradis = np.sqrt(vgrad[0]**2 + vgrad[1]**2)
    if deramp:
        gradis = deramp_unw_np(gradis)  # You should implement the deramp_unw function for NumPy arrays
    return gradis

def deramp_unw_np(array):
    """Deramps unwrapped interferogram represented as a NumPy array."""
    # Assuming _detrend_2d_ufunc can work directly with NumPy arrays.
    # Replace np.nan with 0s for compatibility with the detrending function.
    array_nonan = np.nan_to_num(array)
    
    # Apply the detrending function.
    # This is a simplified call assuming _detrend_2d_ufunc is compatible with NumPy arrays.
    detrended_array = _detrend_2d_ufunc(array_nonan)
    
    return detrended_array

# Example usage:
# array = np.random.rand(100, 100)  # Example NumPy array
# detrended_array = deramp_unw_np(array)
    
def _detrend_2d_ufunc(arr):
    assert arr.ndim == 2, "Input array must be 2D."
    rows, cols = arr.shape
    
    # Prepare the design matrix for a plane: 1, x, and y terms.
    col0 = np.ones(rows * cols)
    col1 = np.repeat(np.arange(rows), cols) + 1
    col2 = np.tile(np.arange(cols), rows) + 1
    G = np.vstack([col0, col1, col2]).T
    
    # Reshape the 2D array into a vector for the observed data.
    d_obs = arr.ravel()
    
    # Use numpy.linalg.lstsq for solving the linear least squares problem.
    # It is more numerically stable than manually inverting matrices.
    m_est, _, _, _ = np.linalg.lstsq(G, d_obs, rcond=None)
    
    # Compute the estimated (fitted) data points and reshape to the original 2D array shape.
    d_est = G.dot(m_est)
    linear_fit = d_est.reshape(arr.shape)
    
    # Return the detrended array.
    return arr - linear_fit


def plot_fwr_bwr(sf_array, fwr, bwr, limit=1, colorbar=0, polys=1, save_path=None):
    '''
    The code is to illustrate the polygons and array in together. It should be improved...
    '''

    
    fig, ax = plt.subplots(figsize=(15, 30))  # Make the plot bigger
    
    # Plot the phase mask data
    cax = ax.imshow(sf_array, cmap='bwr', interpolation='none', aspect='auto')
    # plt.colorbar()
    
    # Highlight FWR polygons
    if polys==1:
        for poly in fwr:
            if poly is not None:
                exterior = poly.exterior
                if exterior:
                    x, y = exterior.xy
                    ax.plot(x, y, color='blue', linewidth=1, alpha=0.7)  # Highlight FWR with blue color
                    # Calculate the centroid for labeling
                    centroid_x = sum(x) / len(x)
                    centroid_y = sum(y) / len(y)
                    if limit ==1:
                        ax.text(centroid_x + 25, centroid_y, 'FWR', ha='center', va='center', fontsize=14, clip_on=True)
    
        # Highlight BWR polygons
        for poly in bwr:
            if poly is not None:
                exterior = poly.exterior
                if exterior:
                    x, y = exterior.xy
                    ax.plot(x, y, color='red', linewidth=1, alpha=0.7)  # Highlight BWR with red color
                    # Calculate the centroid for labeling
                    centroid_x = sum(x) / len(x)
                    centroid_y = sum(y) / len(y)
                    if limit ==1:
                        ax.text(centroid_x+25, centroid_y, 'BWR', ha='center', va='center', fontsize=14, clip_on=True)
        
    # ax.set_title('FWR and BWR Polygons')
    ax.set_xlabel('Range Pixel')
    ax.set_ylabel('Azimuth Pixel')
    if limit == 1:
        ax.set_xlim(2000, 2500)  # Zoom in on the x-axis from 2000 to 2500
        ax.set_ylim(7000, 5000)  # Zoom in on the y-axis from 7000 to 5000 (reversed to keep the origin at the top left)

    if colorbar ==1:
        fig.colorbar(cax, ax=ax, orientation='vertical',fraction=0.046, pad=-0.1)

    if save_path:
        plt.savefig(f'{save_path}.png')
    else:
        plt.show()


def adf_flt(phase_data, kernel_size=5, alpha=0.6, median_kernel_size=3):
    """
    Enhanced adaptive filter inspired by Goldstein's approach for phase images, incorporating a median filter
    on the residuals between the original and initially filtered phase data.
    
    Args:
        phase_data (numpy.ndarray): 2D array of phase data.
        kernel_size (int): Size of the Gaussian kernel for the initial smoothing.
        alpha (float): Factor to control the adaptiveness of the initial filter.
        median_kernel_size (int): Size of the kernel for the median filter applied to residuals.
        
    Returns:
        numpy.ndarray: Final, adjusted filtered phase data.
    """
    # Step 0: change zero to NaN before filtering...
    phase_data = np.where(phase_data == np.nan, 0, phase_data)
    
    # Step 1: Convert phase to complex representation
    complex_data = np.exp(1j * phase_data)

    # Step 2: Generate a Gaussian kernel for initial smoothing
    ax = np.linspace(-(kernel_size - 1) / 2., (kernel_size - 1) / 2., kernel_size)
    gauss = np.exp(-0.5 * np.square(ax) / np.square(kernel_size / alpha))
    kernel = np.outer(gauss, gauss)
    kernel /= kernel.sum()

    # Step 3: Apply initial filtering
    filtered_complex = convolve2d(complex_data, kernel, mode='same', boundary='wrap')
    filtered_phase = np.angle(filtered_complex)

    # Step 4: Calculate and smooth residuals with median filter
    residuals = phase_data - filtered_phase
    smoothed_residuals = medfilt(residuals, median_kernel_size)

    # Step 5: Adjust the filtered phase with smoothed residuals
    adjusted_filtered_phase = filtered_phase + smoothed_residuals

    #preserve the edge
    edge_control=phase_data-adjusted_filtered_phase
    condition = (edge_control == phase_data)
    adjusted_filtered_phase = np.where(condition, phase_data, adjusted_filtered_phase)
    adjusted_filtered_phase = np.where(phase_data == 0, np.nan, adjusted_filtered_phase)
    return adjusted_filtered_phase.astype(np.float32) #np.angle(np.exp(1j * adjusted_filtered_phase))  # Ensure the phase is wrapped properly


def medianfilter_array(arr, ws = 32):
    """use dask median filter on array
    works with both xarray and numpy array
    """
    chunksize = (ws*8, ws*8)
    if type(arr)==type(xr.DataArray()):
        inn = arr.values
    else:
        inn = arr
    arrb = da.from_array(inn, chunks=chunksize)
    arrfilt=ndfilters.median_filter(arrb, size=(ws,ws), mode='reflect').compute()
    if type(arr)==type(xr.DataArray()):
        out = arr.copy()
        out.values = arrfilt
    else:
        out = arrfilt
    return out

###median filtering on residuals, it meanings take time more.
from scipy.ndimage import median_filter
import dask.array as da
import xarray as xr
from dask_image.ndfilters import median_filter
# Make sure you have dask_image installed, or install it using pip install dask_image

def medianfilt_res(arr, ws=96):
    """
    Use Dask median filter on array and apply residual filtering.
    Works with both xarray and numpy array.
    
    Args:
        arr: Input array, either a NumPy array or an xarray DataArray.
        ws (int): Window size for the median filter.
        
    Returns:
        The filtered array with residual adjustment, same type as input.
    """
    print(ws)
    chunksize = (ws*8, ws*8)
    
    # Check if input is xarray DataArray and extract values if so
    is_xarray = isinstance(arr, xr.DataArray)
    if is_xarray:
        inn = arr.values
    else:
        inn = arr
    
    # Convert to Dask array for chunked processing
    arrb = da.from_array(inn, chunks=chunksize)
    
    # Apply median filter using dask_image.ndfilters
    arrfilt = median_filter(arrb, size=(ws, ws), mode='reflect').compute()
    
    # Calculate residuals
    residuals = inn - arrfilt
    
    # Smooth residuals using Dask (note: convert back to Dask array with the same chunksize)
    residuals_dask = da.from_array(residuals, chunks=chunksize)
    smoothed_residuals = median_filter(residuals_dask, size=(ws, ws), mode='reflect').compute()
    
    # Adjust the filtered array with smoothed residuals
    adjusted_arrfilt = arrfilt + smoothed_residuals
    
    # Wrap output in xarray DataArray if input was xarray
    if is_xarray:
        out = arr.copy()
        out.values = adjusted_arrfilt
    else:
        out = adjusted_arrfilt
    
    return out


def median_filter_phase(phase_data, median_kernel_size=7):
    """
    Apply median filtering to phase data.
    
    Args:
        phase_data (numpy.ndarray): 2D array of phase data.
        median_kernel_size (int): Size of the kernel for the median filter.
        
    Returns:
        numpy.ndarray: Median-filtered phase data.
    """
    # Replace NaNs with zeros (if needed)
    #phase_data_nonan = np.where(np.isnan(phase_data), 0, phase_data)
    
    # Apply median filtering
    filtered_phase = medfilt2d(phase_data, kernel_size=median_kernel_size)

    # Apply residual filtering
    residuals = phase_data - filtered_phase
    smoothed_residuals = medfilt2d(residuals, median_kernel_size)

    # Step 5: Adjust the filtered phase with smoothed residuals
    adjusted_filtered_phase = filtered_phase + smoothed_residuals
    
    # Restore NaNs where original data was NaN (if needed)
    #filtered_phase = np.where(np.isnan(phase_data), np.nan, filtered_phase)
    
    return adjusted_filtered_phase.astype(np.float32)


def gaussian_filter_phase(phase_data, sigma=1.0, median_kernel_size=3):
    """
    Apply Gaussian filtering to phase data.
    
    Args:
        phase_data (numpy.ndarray): 2D array of phase data.
        sigma (float): Standard deviation for Gaussian kernel.
        
    Returns:
        numpy.ndarray: Gaussian-filtered phase data.
    """
    # Replace NaNs with zeros (if needed)
    #phase_data_nonan = np.where(np.isnan(phase_data), 0, phase_data)
    
    # Apply Gaussian filtering
    filtered_phase = gaussian_filter(phase_data, sigma=sigma)
    
    # Apply residual filtering
    residuals = phase_data - filtered_phase
    smoothed_residuals = medfilt2d(residuals, median_kernel_size)

    # Step 5: Adjust the filtered phase with smoothed residuals
    adjusted_filtered_phase = filtered_phase + smoothed_residuals
    
    # Restore NaNs where original data was NaN (if needed)
    #filtered_phase = np.where(np.isnan(phase_data), np.nan, filtered_phase)
    
    return adjusted_filtered_phase.astype(np.float32)
  

#################
# def apply_mask_only(data, polygons, reverse_imag=False):
#     """
#     Applies a mask generated for polygons to the data, setting masked areas to np.nan.
#     For a specific condition (e.g., bwr_filtered), reverses the sign of the imaginary part.
    
#     Parameters:
#     - data: The numpy array to be masked. Can be complex.
#     - polygons: A list of shapely Polygon objects used to generate the mask.
#     - reverse_imag: Boolean flag to indicate if the imaginary part's sign should be reversed.
    
#     Returns:
#     - A numpy array with the same type as 'data' where areas covered by the polygons are
#       unchanged, and areas not covered are set to np.nan. If reverse_imag is True, the
#       imaginary part's sign is reversed for the masked areas.
#     """
#     dtype = np.complex64 if np.iscomplexobj(data) else np.float32
#     masked_data = np.full(data.shape, np.nan, dtype=dtype)
    
#     # Generate mask for polygons
#     mask = generate_mask_for_polygons_optimized(data.shape, polygons)
    
#     if reverse_imag:
#         # Apply the mask and reverse the sign of the imaginary part within the masked area
#         imag_part_reversed = np.conj(data[mask])  # Reverse sign of imag part
#         masked_data[mask] = np.real(data[mask]) + 1j * np.imag(imag_part_reversed)
#     else:
#         # Apply the mask without reversing the sign of the imaginary part
#         masked_data[mask] = data[mask]
    
#     return masked_data



def apply_mask_only(data, polygons, reverse_imag=False, phase_factor=-0.53):
    """
    Applies a mask generated for polygons to the data, setting masked areas to np.nan.
    For a specific condition (e.g., bwr_filtered), multiplies the phase by a factor if reverse_imag is True.
    
    Parameters:
    - data: The numpy array to be masked. Can be complex.
    - polygons: A list of shapely Polygon objects used to generate the mask.
    - reverse_imag: Boolean flag to indicate if the phase should be multiplied by the factor.
    - phase_factor: The factor by which to multiply the phase of the data in the masked areas.
    
    Returns:
    - A numpy array with the same type as 'data' where areas covered by the polygons are
      unchanged, and areas not covered are set to np.nan. If reverse_imag is True,
      the phase of the data is multiplied by the given factor in the masked areas.
    """
    dtype = np.complex64 if np.iscomplexobj(data) else np.float32
    masked_data = np.full(data.shape, np.nan, dtype=dtype)
    
    # Generate mask for polygons
    mask = generate_mask_for_polygons_optimized(data.shape, polygons)
    
    if reverse_imag:
        # Extract magnitude and phase for complex data within the mask
        magnitude = np.abs(data[mask])
        phase = np.angle(data[mask])
        
        # Multiply the phase by the given factor
        modified_phase = phase / phase_factor
        
        # Convert back to complex form using the modified phase
        modified_data = magnitude * np.exp(1j * modified_phase)
        
        masked_data[mask] = modified_data
    else:
        # Apply the mask without modifying the phase
        masked_data[mask] = data[mask]
    
    return masked_data

def scaling_before_adf(data, polygons, sf_array):
    """
    This code apply the scaling factor to bwr and fwr phase. 
    data=np.complex64 sboi data
    polygons= bwr or fwr polygon 
    sf_array= scaling factor of sboi
    """
    dtype = np.complex64 if np.iscomplexobj(data) else np.float32
    
    # Generate mask for polygons
    mask = generate_mask_for_polygons_optimized(data.shape, polygons)
    
    dd_pha=np.where(mask, data, np.nan)
    magnitude = np.abs(dd_pha)
    phase=np.angle(dd_pha)
    
    scaled_phase=phase*sf_array
    modified_data = magnitude * np.exp(1j * scaled_phase)
    modified_data = modified_data.astype(dtype)
    
    return modified_data


def pha2cpx(pha, cpx, az_line, width):
    """
    Modifies the real part of a new complex array generated from phase data ('pha')
    by replacing it with the real part of an existing complex array ('cpx').
    
    Parameters:
    - pha: A numpy array containing phase data.
    - cpx: A complex numpy array from which the real part is extracted and reshaped.
    - az_line: The number of azimuth lines for reshaping 'cpx'.
    - width: The width (number of range pixels) for reshaping 'cpx'.
    
    Returns:
    - A complex numpy array generated from 'pha' with its real part replaced by
      the real part from 'cpx', after reshaping 'cpx' to match 'az_line' and 'width'.
    """
    # Reshape 'cpx' to match the expected dimensions and generate a new complex array from 'pha'
    #cpx_reshaped = cpx.reshape(az_line, width)
    cpx_reshaped = cpx.copy()
    new_cpx_from_pha = np.exp(1j * pha)
    new_cpx_from_pha=new_cpx_from_pha.flatten()
    # Ensure 'new_cpx_from_pha' is compatible with the reshaped 'cpx' dimensions
    if new_cpx_from_pha.shape != cpx_reshaped.shape:
        raise ValueError("Shape mismatch between phase-generated complex array and existing complex array")
    
    # Extract the real part of reshaped 'cpx' and the imaginary part of 'new_cpx_from_pha'
    real_part_old_cpx = np.real(cpx_reshaped)
    imaginary_part_new_cpx_from_pha = np.imag(new_cpx_from_pha)
    
    # Combine the extracted real part with the imaginary part of the new complex array
    modified_cpx = real_part_old_cpx + 1j * imaginary_part_new_cpx_from_pha
    modified_cpx=modified_cpx.reshape(az_line, width)
    return modified_cpx



def create_tab_file(epoch, frame_dir, frame, type='RSLC'):
    # Adjust the file name based on the type
    if type == 'RSLC':
        tab_file = os.path.join(frame_dir, 'tab', epoch + 'R_tab')
    elif type == 'SLC':
        tab_file = os.path.join(frame_dir, 'tab', epoch + '_tab')
    else:
        raise ValueError("Unsupported type specified. Choose 'RSLC' or 'SLC'.")

    # Ensure the tab directory exists
    os.makedirs(os.path.dirname(tab_file), exist_ok=True)

    # Check if the tab file does not exist to avoid overwriting
    if not os.path.exists(tab_file):
        if type == 'RSLC':
            inp = os.path.join('RSLC', epoch, epoch)  # Corrected to match the type in the path
            cmd = f"createSLCtab_frame {inp} rslc {frame}"
        elif type == 'SLC':
            inp = os.path.join('SLC', epoch, epoch)  # Corrected to match the type in the path
            cmd = f"createSLCtab_frame {inp} slc {frame}"
        
        # Execute the command and write its output to the tab file
        with open(tab_file, 'w') as file:
            subprocess.run(cmd, shell=True, stdout=file, check=True)

    return tab_file

def framepath_tab(tab_array, frame_directory):
    
    """This function is just the give the full path of the tab file to read correctly from GAMMA."""
    updated_tab_array = []
    for row in tab_array:
        updated_row = [os.path.join(frame_directory, item) for item in row]
        updated_tab_array.append(updated_row)
    return updated_tab_array


def rasterize_polygon_optimized(polygon, shape):
    # Convert polygon points to a Path object
    path = Path(polygon.exterior.coords)

    # Generate a grid of points across the array
    y, x = np.mgrid[:shape[0], :shape[1]]
    points = np.vstack((x.flatten(), y.flatten())).T

    # Use the Path object to test which points are inside the polygon
    mask = path.contains_points(points).reshape(shape)

    return mask

def s1_azfm(r, t0, azp):
  """azfr = s1_azfm(r, t0, azp)
  
  Calculate azimuth FM rate given slant range, reference slant-range delay and the azimuth FM rate polynomial for ScanSAR data
  
  **Arguments:**
  
  * r:    slant range (meters)
  * t0:   reference slant range time for the polynomial (center swath delay in s)
  * azp:  polynomial coefficients
  
  **Output:**
  
  * the function returns the azimuth FM rate"""

  tsr = 2.0 * r / speed_of_light;
  dt = tsr - t0;
  azfr = azp[0] + dt * (azp[1] + dt*(azp[2] + dt*(azp[3] + dt*azp[4])));
  return azfr;

def generate_mask_for_polygons_optimized(shape, polygons):
    # Initialize masks with the same shape as diff_double_mask_pha, filled with False
    mask = np.zeros(shape, dtype=bool)
    y_indices, x_indices = np.indices(shape)

    for polygon in polygons:
        minx, miny, maxx, maxy = polygon.bounds
        # Convert bounds to indices; you might need to adjust this depending on the coordinate system
        minx_idx, miny_idx = int(minx), int(miny)
        maxx_idx, maxy_idx = int(maxx) + 1, int(maxy) + 1
        # Ensure indices are within the array bounds
        minx_idx, miny_idx = max(minx_idx, 0), max(miny_idx, 0)
        maxx_idx, maxy_idx = min(maxx_idx, shape[1]), min(maxy_idx, shape[0])

        # Create a sub-mask for the bounding box area
        sub_mask = np.zeros_like(mask, dtype=bool)
        sub_mask[miny_idx:maxy_idx, minx_idx:maxx_idx] = True

        # Refine the sub-mask by checking points within the bounding box against the polygon
        sub_y, sub_x = np.ogrid[miny_idx:maxy_idx, minx_idx:maxx_idx]
        for x, y in zip(sub_x.flatten(), sub_y.flatten()):
            if polygon.contains(Point(x, y)):
                mask[y, x] = True
            else:
                sub_mask[y, x] = False

        # Combine the sub-mask with the overall mask
        mask |= sub_mask

    return mask

def get_param_gamma(param, parfile, floatt = True, pos = 0):
    a = grep1line(param,parfile).split()[1+pos]
    if floatt:
        a = float(a)
    return a


def get_dfDC(path_to_slcdir, f0=5405000500, burst_interval = 2.758277, returnka = True, returnperswath = False, returnscalefactor=True):
    #f0 = get_param_gamma('radar_frequency', parfile)
    #burst_interval = get_param_gamma('burst_interval', topsparfile)
    epoch = os.path.basename(path_to_slcdir)
    print(epoch)
    frame = path_to_slcdir.split('/')[-3]
    
    if len(frame)!=17:
        frame=path_to_slcdir.split('/')[-4]
    parfile = os.path.join(path_to_slcdir, epoch+'.rslc.par')
    #parfile = glob.glob(path_to_slcdir+'/????????.slc.par')[0]
    #topsparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.TOPS_par')
    #iwparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.par')
    #
    
    lam = speed_of_light / f0
    dfDC = []
    kas = []
    ctr_range = []
    far_range = []
    near_range = []
    scalefactor= []
    afmrate_srdly = []
    afmrate_ply= []
    kr_list = []
    numbursts = [ int(frame.split('_')[2][:2]), int(frame.split('_')[2][2:4]), int(frame.split('_')[2][4:6])]
    azps_list = []
    az_line_time_list = []
    #krs = []
    #print('This is a proper solution but applied to primary SLC image. originally it is applied by GAMMA on the RSLC...')
    #for n in range(len(topsparfiles)):
    for n in [1,2,3]:
        topsparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n)+'.rslc.TOPS_par')
        iwparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n)+'.rslc.par')
        if (not os.path.exists(iwparfile)) or (not os.path.exists(topsparfile)):
            dfDC.append(np.nan)
            kas.append(np.nan)
            numbursts[n-1] = np.nan
        else:
            #topsparfile = topsparfiles[n]
            #iwparfile = iwparfiles[n]
            az_steering_rate = get_param_gamma('az_steering_rate', topsparfile) # az_steering_rate is the antenna beam steering rate
            #az_ster_rate.append(az_steering_rate)
            r1 = get_param_gamma('center_range_slc', iwparfile)
            #get the satellite velocity
            midNstate = int(get_param_gamma('number_of_state_vectors', iwparfile)/2)+1
            sv = 'state_vector_velocity_' + str(midNstate)
            velc1 = get_param_gamma(sv, iwparfile, pos=0)
            velc2 = get_param_gamma(sv, iwparfile, pos=1)
            velc3 = get_param_gamma(sv, iwparfile, pos=2)
            vsat = np.sqrt(velc1**2 + velc2**2 + velc3**2)
            midNstate=1
            # now some calculations
            afmrate_srdelay = get_param_gamma('az_fmrate_srdelay_'+ str(midNstate), topsparfile)
            afmrate_srdly.append(afmrate_srdelay) 
            afmrate_poly = []
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 0))
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 1))
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 2))
            try:
                afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 3))
            except:
                afmrate_poly.append(0)
            try:
                afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 4))
            except:
                afmrate_poly.append(0)
            afmrate_ply.append(afmrate_poly)
            ka = s1_azfm(r1, afmrate_srdelay, afmrate_poly) #unit: Hz/s == 1/s^2
            kr = -2.0 * vsat * az_steering_rate*(pi / 180.0) / lam
            kr_list.append(kr)
            if (kr != 0.0):
                #kt = ka * kr/(kr - ka)
                # but maybe should be kt = (kr*ka)/(ka-kr) # see https://iopscience.iop.org/article/10.1088/1755-1315/57/1/012019/pdf  --- and remotesensing-12-01189-v2, and Fattahi et al...
                # ok, gamma reads kr to be ka... corrected
                kt = kr * ka/(ka - kr)
            else:
                kt = -ka
            #finally calculate dfDC:
            #burst_interval = get_param_gamma('burst_interval', topsparfile)
            kas.append(ka)
            #krs.append(kr)
            dfDC.append(kt*burst_interval) #burst_interval is time within the burst... we can also just calculate.. see Grandin: eq 15: hal.archives-ouvertes.fr/hal-01621519/document
            #ok, that's the thing - burst_interval is actually t(n+1) - t(n) - see remotesensing-12-01189-v2
            #so it should be kt * -burst_interval, that is why GAMMA has the -kt J ... ok, good to realise this

            ####calculating scalefactor (Nergizci)
            azps=np.float64(get_param_gamma('azimuth_pixel_spacing', iwparfile))
            if not azps_list:
                azps_list.append(azps)
            az_line_time=np.float64(get_param_gamma('azimuth_line_time', iwparfile))
            if not az_line_time_list:
                az_line_time_list.append(az_line_time)
            
            dfdc=kt*burst_interval
            sf=(azps)/(dfdc*az_line_time*2*np.pi)
            scalefactor.append(sf)
            
            ###calculating ssd (Nergizci)
            ctr_range_temp=get_param_gamma('center_range_slc', iwparfile)
            far_range_temp=get_param_gamma('far_range_slc', iwparfile)
            near_range_temp=get_param_gamma('near_range_slc', iwparfile)
            ctr_range.append(ctr_range_temp)
            far_range.append(far_range_temp)
            near_range.append(near_range_temp)
    
    #print(scalefactor)
    #for sw1 and sw2
    
    r_so=(far_range[0] + near_range[1])/2
    ka_sb_list=[]
    kt_sb_list=[]
    t_bd_list=[]
    for n in range(2):
        topsparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n+1)+'.rslc.TOPS_par')
        iwparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n+1)+'.rslc.par')
        ka_sb=s1_azfm(r_so, afmrate_srdly[n], afmrate_ply[n])
        ka_sb_list.append(ka_sb)
        if (kr_list[n] != 0.0):
            kt_sb = ka_sb_list[n]* kr_list[n]/(ka_sb_list[n] - kr_list[n])
        else:
            kt_sb = -ka_sb_list[n]
        kt_sb_list.append(kt_sb)

        t_bd_list_temp = []
        #get burst delta times
        for i in range(numbursts[n]):
            t_bd=get_param_gamma(f'burst_start_time_{i+1}', topsparfile)
            t_bd_list_temp.append(t_bd)

        t_bd_list.append(t_bd_list_temp)

    # Directly calculate dt_f for each corresponding pair across the two lists
    dt_f_list = [t_bd_list[1][i] - t_bd_list[0][i] for i in range(len(t_bd_list[0]))]
    
    dt_b_list = []
    for dt_f in dt_f_list:
        if dt_f < 0:
            dt_b = burst_interval + dt_f
        else:
            dt_b = dt_f - burst_interval
        dt_b_list.append(dt_b)
    
    # Now calculate dfDC_f and dfDC_b using the time differences
    dfDC_f_list = [(kt_sb_list[0]+kt_sb_list[1])*dt_f/2.0 for dt_f in dt_f_list]
    dfDC_b_list = [(kt_sb_list[0]+kt_sb_list[1])*dt_b/2.0 for dt_b in dt_b_list]

    # ####calculating scalefactor
    
    if len(azps_list) > 0 and len(az_line_time_list) > 0:
        azps = azps_list[0]
        az_line_time = az_line_time_list[0]

        # Calculate sff and sfb for each dfDC_f and dfDC_b value
        sff1 = [azps / (df * -1*az_line_time * 2 * np.pi) for df in dfDC_f_list]
        sfb1 = [azps / (db * -1*az_line_time * 2 * np.pi) for db in dfDC_b_list]
    
    else:
        print("azps_list or az_line_time_list is empty, cannot calculate sff and sfb.")

    ####################################################
    #for subswath 2-3
    r_so=(far_range[1] + near_range[2])/2
    ka_sb_list=[]
    kt_sb_list=[]
    t_bd_list=[]
    for n in range(1,3,1):
        topsparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n+1)+'.rslc.TOPS_par')
        iwparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n+1)+'.rslc.par')
        ka_sb=s1_azfm(r_so, afmrate_srdly[n], afmrate_ply[n])
        ka_sb_list.append(ka_sb)
        if (kr_list[n] != 0.0):
            kt_sb = ka_sb_list[n-1]* kr_list[n]/(ka_sb_list[n-1] - kr_list[n])
        else:
            kt_sb = -ka_sb_list[n-1]
        kt_sb_list.append(kt_sb)

        t_bd_list_temp = []
        #get burst delta times
        for i in range(numbursts[n]):
            t_bd=get_param_gamma(f'burst_start_time_{i+1}', topsparfile)
            t_bd_list_temp.append(t_bd)

        t_bd_list.append(t_bd_list_temp)


    # Directly calculate dt_f for each corresponding pair across the two lists
    dt_f_list = [t_bd_list[1][i] - t_bd_list[0][i] for i in range(len(t_bd_list[0]))]

    dt_b_list = []
    for dt_f in dt_f_list:
        if dt_f < 0:
            dt_b = burst_interval + dt_f
        else:
            dt_b = dt_f - burst_interval
        dt_b_list.append(dt_b)
    
    # Now calculate dfDC_f and dfDC_b using the time differences
    dfDC_f_list = [(kt_sb_list[0]+kt_sb_list[1])*dt_f/2.0 for dt_f in dt_f_list]
    dfDC_b_list = [(kt_sb_list[0]+kt_sb_list[1])*dt_b/2.0 for dt_b in dt_b_list]

    # ####calculating scalefactor
    
    if len(azps_list) > 0 and len(az_line_time_list) > 0:
        azps = azps_list[0]
        az_line_time = az_line_time_list[0]

        # Calculate sff and sfb for each dfDC_f and dfDC_b value
        sff2 = [azps / (df *1* az_line_time * 2 * np.pi) for df in dfDC_f_list]
        sfb2 = [azps / (db *1*az_line_time * 2 * np.pi) for db in dfDC_b_list]
    
    else:
        print("azps_list or az_line_time_list is empty, cannot calculate sff and sfb.")
    
    # Create a list with two columns, where each row is [sff1[i], sff2[i]]
    final_sff_temp = [[sff1[i], sff2[i]] for i in range(len(sff1))]
    # Reshape sff from 25x2 to a 50x1 list
    final_sff = [row[0] for row in final_sff_temp] + [row[1] for row in final_sff_temp]
    
    final_sfb_temp = [[sfb1[i], sfb2[i]] for i in range(len(sfb1))]
    # Reshape sff from 25x2 to a 50x1 list
    final_sfb = [row[0] for row in final_sfb_temp[1:]] + [row[1] for row in final_sfb_temp[1:]]
    if not returnperswath:
        numbursts = np.array(numbursts)
        dfDC = np.nansum(numbursts*np.array(dfDC)) / np.sum(numbursts)
        ka = np.nansum(numbursts*np.array(kas)) / np.sum(numbursts)
    #kr = np.mean(krs)
    if returnka:
        return dfDC, ka #, kr
    if returnscalefactor:
        return scalefactor, final_sff, final_sfb
    else:
        return dfDC

##### scaling factor in pixel base.
def get_sf_array(path_to_slcdir, f0=5405000500, burst_interval = 2.758277):
    #f0 = get_param_gamma('radar_frequency', parfile)
    #burst_interval = get_param_gamma('burst_interval', topsparfile)
    epoch = os.path.basename(path_to_slcdir)
    frame = path_to_slcdir.split('/')[-3]
    parfile = os.path.join(path_to_slcdir, epoch+'.rslc.par')
    #parfile = glob.glob(path_to_slcdir+'/????????.slc.par')[0]
    #topsparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.TOPS_par')
    #iwparfiles = glob.glob(path_to_slcdir+'/????????.IW?.slc.par')
    #
    
    lam = speed_of_light / f0
    numbursts = [ int(frame.split('_')[2][:2]), int(frame.split('_')[2][2:4]), int(frame.split('_')[2][4:6])]
    # print(numbursts)
    azps_list = []
    az_line_time_list = []
    rps_list=[]
    ka_dict={}
    afmrate_srdelay_dict = {}
    afmrate_poly_dict= {}
    ext_burst_win_r0_dict={}
    ext_burst_win_r1_dict={}
    ext_burst_win_a0_dict={}
    ext_burst_win_a1_dict={}
    burst_start_time_dict={}
    az_line_time_dict={}
    lpb_dict={}
    far_range_dict={}
    near_range_dict={}
    ctr_range_dict={}
    #print('This is a proper solution but applied to primary SLC image. originally it is applied by GAMMA on the RSLC...')
    for n in [1,2,3]:
        key=f'SW{n}'
        topsparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n)+'.rslc.TOPS_par')
        iwparfile = os.path.join(path_to_slcdir, epoch+'.IW'+str(n)+'.rslc.par')
        # print(topsparfile, iwparfile)
        if (not os.path.exists(iwparfile)) or (not os.path.exists(topsparfile)):
            dfDC.append(np.nan)
            kas.append(np.nan)
            numbursts[n-1] = np.nan
        else:
            #topsparfile = topsparfiles[n]
            #iwparfile = iwparfiles[n]
            az_steering_rate = get_param_gamma('az_steering_rate', topsparfile) # az_steering_rate is the antenna beam steering rate
            #az_ster_rate.append(az_steering_rate)
            r1 = get_param_gamma('center_range_slc', iwparfile)
            #get the satellite velocity
            midNstate = int(get_param_gamma('number_of_state_vectors', iwparfile)/2)+1
            sv = 'state_vector_velocity_' + str(midNstate)
            velc1 = get_param_gamma(sv, iwparfile, pos=0)
            velc2 = get_param_gamma(sv, iwparfile, pos=1)
            velc3 = get_param_gamma(sv, iwparfile, pos=2)
            vsat = np.sqrt(velc1**2 + velc2**2 + velc3**2)
            midNstate=1
            # now some calculations
            afmrate_srdelay = get_param_gamma('az_fmrate_srdelay_'+ str(midNstate), topsparfile)
            afmrate_srdelay_dict[key]=afmrate_srdelay 
            afmrate_poly = []
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 0))
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 1))
            afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 2))
            try:
                afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 3))
            except:
                afmrate_poly.append(0)
            try:
                afmrate_poly.append(get_param_gamma('az_fmrate_polynomial_' + str(midNstate), topsparfile, pos = 4))
            except:
                afmrate_poly.append(0)
            afmrate_poly_dict[key]=afmrate_poly

            ##doppler rate calculation
            kr = -2.0 * vsat * az_steering_rate*(pi / 180.0) / lam  # ka for GAMMA
            ka_dict[key]=kr 
            
  
            for i in range(numbursts[n-1]):
                keys = f'SW{n}_{i+1}' ##key for dictionary?!
                ext_burst_win_r0=get_param_gamma('burst_win_' +str(i+1), topsparfile, pos=0)
                ext_burst_win_r0_dict[keys] = ext_burst_win_r0
                ext_burst_win_r1=get_param_gamma('burst_win_' +str(i+1), topsparfile, pos=1)
                ext_burst_win_r1_dict[keys] = ext_burst_win_r1
                ext_burst_win_a0=get_param_gamma('burst_win_' +str(i+1), topsparfile, pos=2)
                ext_burst_win_a0_dict[keys]=ext_burst_win_a0
                ext_burst_win_a1=get_param_gamma('burst_win_' +str(i+1), topsparfile, pos=3)
                ext_burst_win_a1_dict[keys]=ext_burst_win_a1
                burst_start_time=get_param_gamma('burst_start_time_'+str(i+1), topsparfile, pos=0)
                burst_start_time_dict[keys]=burst_start_time
            
            
            
            azps=np.float64(get_param_gamma('azimuth_pixel_spacing', iwparfile))
            rps=np.float64(get_param_gamma('range_pixel_spacing', iwparfile))
            rps_ml=rps*20 ##multilooked pixel of range
            azps_list.append(azps)
            rps_list.append(rps)
            
            az_line_time=np.float64(get_param_gamma('azimuth_line_time', iwparfile))
            az_line_time_dict[key]=az_line_time

            
                
            ctr_range_temp=get_param_gamma('center_range_slc', iwparfile)
            ctr_range_dict[key]=ctr_range_temp
            far_range_temp=get_param_gamma('far_range_slc', iwparfile)
            far_range_dict[key]=far_range_temp
            near_range_temp=get_param_gamma('near_range_slc', iwparfile)
            near_range_dict[key]=near_range_temp

            lpb=get_param_gamma('lines_per_burst', topsparfile)
            lpb_dict[key]=lpb

    #print(afmrate_poly_dict)
    #print(afmrate_srdelay_dict)
    # print(ka_dict)
    #print(ext_burst_win_r0_dict)
    #seems everything is fine till here! 
    
    #######regarding the set maxima on subswath-2 but it would be more efficent?
    ##### find about the variation in range for subswath overlaps.
    r_so_results={}
    x_coord={}

############ Set the subswath where the maximum burst exist.
    # Step 1: Determine the sub-swath with the maximum burst number
    burst_count = {}
    
    for key in ext_burst_win_r1_dict:
        parts = key.split('_')
        sw_number = int(parts[0][-1])
        burst_number = int(parts[-1])
        
        if sw_number not in burst_count:
            burst_count[sw_number] = burst_number
        else:
            if burst_number > burst_count[sw_number]:
                burst_count[sw_number] = burst_number
    
    # Find the sub-swath with the maximum burst number
    max_sw_number = max(burst_count, key=burst_count.get)

    # Step 2: Collect the burst numbers for the sub-swath with the maximum burst number
    burst_numbers_sw_max = []
    
    for key in ext_burst_win_r1_dict:
        parts = key.split('_')
        sw_number = int(parts[0][-1])
        burst_number = int(parts[-1])
        
        if sw_number == max_sw_number:
            burst_numbers_sw_max.append(burst_number)
    
    # Ensure burst numbers are unique and sorted for proper matching
    burst_numbers_sw_max = sorted(set(burst_numbers_sw_max))

    # print(ext_burst_win_a0_dict)
    #print(ext_burst_win_r1_dict)


    for key, value in ext_burst_win_r1_dict.items():
        parts = key.split('_')
        sw_number= int(parts[0][-1]) # Correctly extract sw number and burst number
        # Assume keys in ext_burst_win_r0_dic are in the same format and case-sensitive
        if sw_number == 1 or sw_number == 2:  # For subswath1's and subswath2's in r1_dic
            for burst_number in burst_numbers_sw_max:
                next_sw_key_r0 = f'SW{sw_number+1}_{burst_number}'  # Correct key for matching in r0_dic
                # print(next_sw_key_r0)
                if next_sw_key_r0 in ext_burst_win_r0_dict:
                    # print(f'inside the ext {next_sw_key_r0}')
                    if int(next_sw_key_r0.split('_')[0][2]) == int(key.split('_')[0][2]) + 1 and int(next_sw_key_r0.split('_')[1]) == int(key.split('_')[1]):    ## We have taken the sw2 burst number setting and then make elemination here to catch the matched burst in SW1-2 and SW2-3, valid.
                        r0 = ext_burst_win_r0_dict[next_sw_key_r0]  # r0 is SW2_burst for SW1-2 overlap and SW3_burt for SW2-3, make sense!
                        r0_ml=((r0-rps/2)+rps_ml/2)-250 ##optimize the multilooked center location of pixel 
                        r1= ext_burst_win_r1_dict[key] # or r1 = value r1 is SW1_burst for SW1-2 overlap and SW2_burst for SW2-3, another sensible moment from MN! 
                        r1_ml=((r1-rps/2)+rps_ml/2)+250  ##same reason;multilooked center pixel. 
                        nr=int(np.ceil(r1_ml-r0_ml)/rps_ml)+1 ##number of pixel ib
                        r_so = r0_ml + np.arange(nr) * rps_ml
                        r_so_results[next_sw_key_r0] = r_so ###recheck the calculation of r_so_result!
                        x_coord[next_sw_key_r0]=r_so
                        last_r_so = r_so
                    else:
                        # If burst number not coincide with the subswaths? It is not make sense to me valid this for each case, but works fine so far! investigate regarding feedback?
                        if last_r_so is not None:
                            r_so_results[next_sw_key_r0] = last_r_so
                            x_coord[next_sw_key_r0] = last_r_so
                else:
                    # print(f'outside {next_sw_key_r0}')
                    # If next_sw_key_r0 is not in ext_burst_win_r0_dict, use the last valid r_so
                    if last_r_so is not None:
                        r_so_results[next_sw_key_r0] = last_r_so
                        x_coord[next_sw_key_r0] = last_r_so
    
    def calculate_kr_and_kt(r_so, afmrate_srdelay, afmrate_poly, ka):
        """
        Calculate kr and kt values based on given parameters as pixel based in range direction.
    
        Args:
            r_so (array): Range sequence over which to calculate the parameters.
            afmrate_srdelay (float): Azimuth FM rate slant-range delay.
            afmrate_poly (list): Coefficients of the azimuth FM rate polynomial.
            ka (float): Doppler frequency rate due to satellite motion.
    
        Returns:
            tuple: The calculated kr and kt values.
        """
        kr = s1_azfm(r_so, afmrate_srdelay, afmrate_poly)
        if ka != 0.0:
            kt = ka * kr / (kr - ka)
        else:
            kt = -kr
        return kr, kt
        
    def process_r_so_results(r_so_results, sw_numbers, afmrate_srdelay_dict, afmrate_poly_dict, ka_dict):
        """
        loop calculation of kr and kt values for each burst using calculate_kr_and_kt.
    
        Args:
            r_so_results (dict): Dictionary containing range subswath overlaps.
            sw_numbers (list): List of subswath numbers being processed [1, 2] or [2,3].
            afmrate_srdelay_dict (dict): Dictionary of azimuth FM rate slant-range delays for each subswath.
            afmrate_poly_dict (dict): Dictionary of azimuth FM rate polynomial coefficients for each subswath.
            ka_dict (dict): Dictionary of Doppler frequency rates due to satellite motion for each subswath.
    
        Returns:
            dict, dict: Dictionaries of kr and kt values for each burst within the subswath overlaps.
        """
        kr_dict, kt_dict = {}, {}
        for key, r_so in r_so_results.items():
            if key.startswith(f"SW{sw_numbers[-1]}"):
                #print(key)
                burst_id=key.split('_')[-1]
                for sw_number in sw_numbers:
                    afmrate_srdelay = afmrate_srdelay_dict[f'SW{sw_number}']
                    afmrate_poly = afmrate_poly_dict[f'SW{sw_number}']
                    ka = ka_dict[f'SW{sw_number}']
                    kr, kt = calculate_kr_and_kt(r_so, afmrate_srdelay, afmrate_poly, ka)
                    kr_key = f"SW{sw_number}_{burst_id}"
                    kr_dict[kr_key] = kr
                    kt_dict[kr_key] = kt
        return kr_dict, kt_dict

    def calculate_azi_line(path_to_slcdir, ext_burst_win_a0_dict, ext_burst_win_a1_dict, sw_numbers, azps, kr_dict, kt_dict):
        y_coord_dict={}
        sf_dict={}
        parfile = os.path.join(path_to_slcdir, epoch+'.rslc.par')
        par = pg.ParFile(parfile)
        for i in range(len(sw_numbers)-1):
            sw1, sw2 = sw_numbers[i], sw_numbers[i+1]
            ###azimuth multilooking!
            tazi=par.get_value('azimuth_line_time', index = 0, dtype = float)
            tazi_ml=tazi*4 ##multilooking 4 in azimuth
            
             # Extract keys for both subswaths. Just to give the keys for each overlap. It callced sw1-sw2 but it is cover sw2-sw3 in same time because of loop.
            sw1_keys_full = [key for key in burst_start_time_dict.keys() if key.startswith(f"SW{sw1}_")]
            sw2_keys_full = [key for key in burst_start_time_dict.keys() if key.startswith(f"SW{sw2}_")]

            ###here is totally wrong in case the subswath is not start with same azimuth line time! 
            ###the idea should be find the shift between SW1-2 and SW2-3 if there is using ext_burst_win_a0!
            min_length = min(len(sw1_keys_full), len(sw2_keys_full)) ##it looks logical to me tbh. You can overlap at most the minumum burst number between SW1-2 and SW2-3. Keep going
            # Trim both lists to the minimum length
            sw1_keys = sw1_keys_full[:min_length]
            sw2_keys = sw2_keys_full[:min_length]

            ###Try to find the offset between subswath for irregular pattern, like Sulawesi and Tien-Shan extreme frame: '134D_08966_081302' and '056A_04947_282019'
            for key, value in ext_burst_win_a0_dict.items():
                if key == f'SW{sw1}_1':
                    sw1_a0 = value
                if key == f'SW{sw2}_1':
                    sw2_a0 = value

            ####
            if sw1 == '1':
                offset_a0 = round((sw2_a0 - sw1_a0) / 3.7) ##azimuth line time is around 2.75 between adjacent vertical burst therefore mod3.7 (3.7 because I increase 1 offset could be 1.8 or 0.9 even they match well)gives you how many burst offset exist between start burst in different subswath.?! 
                if offset_a0 > 0:
                    sw1_keys = sw1_keys_full[offset_a0:]  # Adjust keys for sw1 based on the offset
                    sw2_keys = sw2_keys_full[:]           # Keep keys for sw2 as they are
                elif offset_a0 < 0:
                    sw1_keys = sw1_keys_full[:]
                    sw2_keys = sw2_keys_full[np.abs(offset_a0):]
                else:
                    sw1_keys = sw1_keys_full[:]
                    sw2_keys = sw2_keys_full[:]   

            if sw1 == '2':
                offset_a0 = round((sw2_a0 - sw1_a0) / 3.7)
                if offset_a0 > 0:
                    sw1_keys = sw1_keys_full[offset_a0:]  # Adjust keys for sw1 based on the offset
                    sw2_keys = sw2_keys_full[:]           # Keep keys for sw2 as they are
                elif offset_a0 < 0:
                    sw1_keys = sw1_keys_full[:]
                    sw2_keys = sw2_keys_full[np.abs(offset_a0):]
                else:
                    sw1_keys = sw1_keys_full[:]
                    sw2_keys = sw2_keys_full[:]  

            for sw1_key in sw1_keys:
                for sw2_key in sw2_keys:
                    ###interpolation incase kt[sw_1] and kt[sw_2] is not same.
                    if kt_dict[sw1_key].shape[0] != kt_dict[sw2_key].shape[0]:
                        # Determine which array is shorter
                        if kt_dict[sw1_key].shape[0] < kt_dict[sw2_key].shape[0]:
                            shorter_key, longer_key = sw1_key, sw2_key
                        else:
                            shorter_key, longer_key = sw2_key, sw1_key
                            
                        # Original index range for the shorter array
                        x_original = np.arange(kt_dict[shorter_key].shape[0])
                        
                        # New index range to match the length of the longer kt array
                        x_new = np.linspace(0, x_original[-1], kt_dict[longer_key].shape[0])  
                        # Create interpolation function for the shorter array
                        f_interp = interp1d(x_original, kt_dict[shorter_key], kind='linear', fill_value="extrapolate")

                        # Apply interpolation to match the length
                        kt_dict[shorter_key] = f_interp(x_new)
                        
                    if ext_burst_win_a1_dict[sw1_key]>ext_burst_win_a0_dict[sw2_key] and ext_burst_win_a1_dict[sw2_key]>ext_burst_win_a0_dict[sw1_key]: ### check burst overlaps: bwr and fwr!
                        # print(sw1_key, sw2_key)
                        t0_sw1=burst_start_time_dict[sw1_key]
                        t1_sw1=t0_sw1+az_line_time_dict[f'SW{sw1}']*(lpb_dict[f'SW{sw1}']-1)/2
                        t0_sw2=burst_start_time_dict[sw2_key]
                        t1_sw2=t0_sw2+az_line_time_dict[f'SW{sw2}']*(lpb_dict[f'SW{sw2}']-1)/2
                        #### calculation of overlap interval
                        if ext_burst_win_a0_dict[sw2_key] > ext_burst_win_a0_dict[sw1_key]: ### situation of forward SOI
                            t0=ext_burst_win_a0_dict[sw2_key]
                            t0_ml=(t0-tazi/2)+tazi_ml/2    ##Applied same pixel_center convention for azi?
                            t2=ext_burst_win_a1_dict[sw1_key]
                            t2_ml=(t2-tazi/2)+tazi_ml/2
                            
                            naz = int(np.ceil((t2_ml - t0_ml)/tazi_ml))
                            ti_ml = t0_ml + np.arange(naz) * tazi_ml
                            y_coord_dict[f'fwr{sw2_key}']=ti_ml
                            # vector multiplication to form a matrix with range and azimuth dimensions
                            t1_sw1_ml=(t1_sw1-tazi/2)+tazi_ml/2
                            t1_sw2_ml=(t1_sw2-tazi/2)+tazi_ml/2
                            
                            fDC_sw1 = np.outer(kt_dict[sw1_key], (ti_ml - t1_sw1_ml))
                            fDC_sw2 = np.outer(kt_dict[sw2_key], (ti_ml - t1_sw2_ml))
                            dfDC = fDC_sw1 - fDC_sw2
                            sf = (azps) / (dfDC * tazi *2* np.pi)
                            sf_dict[f'fwr{sw2_key}']=sf
                        else:
                            t0=ext_burst_win_a0_dict[sw1_key]
                            t0_ml=(t0-tazi/2)+tazi_ml/2 
                            t2=ext_burst_win_a1_dict[sw2_key]
                            t2_ml=(t2-tazi/2)+tazi_ml/2
                            
                            naz = int(np.ceil((t2_ml - t0_ml)/tazi_ml))
                            ti_ml = t0_ml + np.arange(naz) * tazi_ml
                            y_coord_dict[f'bwr{sw2_key}']=ti_ml
                            t1_sw1_ml=(t1_sw1-tazi/2)+tazi_ml/2
                            t1_sw2_ml=(t1_sw2-tazi/2)+tazi_ml/2
    
                            fDC_sw1 = np.outer(kt_dict[sw1_key], (ti_ml - t1_sw1_ml))
                            fDC_sw2 = np.outer(kt_dict[sw2_key], (ti_ml - t1_sw2_ml))
                            
                            dfDC = fDC_sw1 - fDC_sw2
                            sf = (azps) / (dfDC * tazi *2* np.pi)
                            sf_dict[f'bwr{sw2_key}']=sf
            return y_coord_dict, sf_dict

    # Initialize sw_numbers for SW1 and SW2, then update it for SW2 and SW3
    sw_numbers = ['1', '2']
    # Check if '1' is in sw_numbers and perform calculations
    if '1' in sw_numbers:
        kr_dict, kt_dict = process_r_so_results(r_so_results, sw_numbers, afmrate_srdelay_dict, afmrate_poly_dict, ka_dict)
        y_coord1, sf_dict1 = calculate_azi_line(path_to_slcdir, ext_burst_win_a0_dict, ext_burst_win_a1_dict, sw_numbers, azps, kr_dict, kt_dict)
    else:
        y_coord1, sf_dict1 = np.nan, {}
    
    # Update sw_numbers for SW2 and SW3
    sw_numbers = ['2', '3']
    # Check if '3' is in sw_numbers and perform calculations
    if '3' in sw_numbers:
        kr_dict, kt_dict = process_r_so_results(r_so_results, sw_numbers, afmrate_srdelay_dict, afmrate_poly_dict, ka_dict)
        y_coord2, sf_dict2 = calculate_azi_line(path_to_slcdir, ext_burst_win_a0_dict, ext_burst_win_a1_dict, sw_numbers, azps, kr_dict, kt_dict)
    else:
        y_coord2, sf_dict2 = np.nan, {}

    ###reverse sw1 values:
    sf_dict1= {key: value * -1 for key, value in sf_dict1.items()}
    ##merging SW1-SW2 SOI and SW2-SW3 SOI
    sf_dict={**sf_dict1, **sf_dict2}
    y_coord={**y_coord1, **y_coord2}
    
    
    ###define the path:
    prime=os.path.basename(path_to_slcdir)
    RSLC_dir=path_to_slcdir
    
    # number of subswaths
    nsubswaths = 3
    
    # number of looks
    rlks = 20
    azlks = 4
    # parameter file of the mosaic
    SLC_par = pg.ParFile(os.path.join(RSLC_dir, prime + '.rslc.par'))
    
    # # read SLC parameters
    r0 = SLC_par.get_value('near_range_slc', index = 0, dtype = float)        # near range distance
    rps = SLC_par.get_value('range_pixel_spacing', index = 0, dtype = float)  # range pixel spacing
    t0 = SLC_par.get_value('start_time', index = 0, dtype = float)            # time of first azimuth line
    tazi = SLC_par.get_value('azimuth_line_time', index = 0, dtype = float)   # time between each azimuth line
   
    y_coord_pix={}
    for key, y_values in y_coord.items():
        y_coord_pix[key]= [(round((y - t0)/tazi))/azlks for y in y_values]
    x_coord_pix={}
    for key, x_values in x_coord.items():
        x_coord_pix[key]= [(round((x - r0)/rps))/rlks for x in x_values]


    ##Reading Phase data and its shape to define our matrix
    mli_pardir=pg.ParFile(os.path.join(RSLC_dir, prime +'.rslc.mli.par'))
    width=mli_pardir.get_value('range_samples', index = 0, dtype = int)
    length=mli_pardir.get_value('azimuth_lines', index=0, dtype = int)
    ##defining dummy array to fulfill the values with sf thanks to x_coordinates and y_coordinates
    array_dummy=np.zeros((length, width))

    for key in sf_dict.keys():
        x_key=key[3:]

        if x_key in x_coord_pix:
            x_pixels=np.array(x_coord_pix[x_key])
            y_pixels=np.array(y_coord_pix[key])
            
            sf_values=np.array(sf_dict[key])
            # print(np.shape(sf_values))
            # print(np.shape(x_pixels))
            # print(np.shape(y_pixels))
            for xi in range(len(x_pixels)):
                for yi in range(len(y_pixels)):
                    
                    ##Extract the current scale value
                    sf_value=sf_values[xi, yi]

                    ##calculate the x and y indices in the array dummy
                    x_index=int(x_pixels[xi])
                    y_index=int(y_pixels[yi])


                    # Check if the indices are within the bounds of the array
                    if 0 <= x_index < width and 0 <= y_index < length:
                        # Assign the scale factor value to the correct position in array_dummy
                        array_dummy[y_index, x_index] = sf_value
                        #sf_array=array_dummy.copy()
    # print(np.shape(array_dummy))
    return array_dummy