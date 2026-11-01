#!/usr/bin/env python
"""
Function to extract NISAR frame boundary boxes and save them as polygon files
to be used by the LiCSAR portal R scripts.
"""

import os
import numpy as np
import geopandas as gpd
from shapely.geometry import box, mapping
from rasterio.plot import show
from rasterio.warp import transform_bounds
import xarray as xr
import h5py
from LiCSAR_lib.nisardata import load_gunw


def get_nisar_frame_from_file(h5file):
    """
    Extract frame ID from NISAR file name
    e.g., NISAR_L2_PR_GUNW_009_034_A_018_4005_DHDH_A_20251230T130752_20251230T130827_X05009_N_F_J_001.h5
    Returns: (band, track, path, frame, polarization)
    """
    parts = os.path.basename(h5file).split('_')
    
    # Parse filename to get frame info
    band = parts[1][0]  # 'L' or 'S' - L for GUNW, S for GSLC
    if band == 'G':  # GUNW files are named differently
        band = parts[1]
    else:
        band = parts[2][0]  # for GSLC files
    
    # Determine product type (GUNW or GSLC)
    prodtype = parts[3]
    
    # Parse frame identification
    # Frame pattern: A/D (ascending/descending), track (3-digit), path (3-digit), frame (6-digit)
    frame_part = parts[8]  # e.g., 4005_DHDH_A
    trk = int(frame_part[8:11])  # track number
    path = frame_part[11:14]  # path number
    frame = frame_part[14:20]  # frame number (6 digits)
    polar = frame_part[20]  # polarization (A/D)
    
    return band, trk, path, frame, polar


def extract_frame_bbox_from_gunw(h5file, output_geojson=True):
    """
    Extract the boundary box (bounding polygon) from a NISAR GUNW file
    
    Args:
        h5file: Path to GUNW NetCDF file
        output_geojson: If True, output GeoJSON; if False, output WKT
    
    Returns:
        shapely Polygon of the frame boundaries
    """
    if not os.path.exists(h5file):
        print(f"Error: File {h5file} does not exist")
        return None
    
    # Load the GUNW data to get the geospatial extent
    try:
        # Load GUNW data - we only need the spatial coordinates and CRS
        xrds = xr.open_dataset(h5file)
        
        # Get coordinate arrays
        x_coords = xrds.x.values if 'x' in xrds else xrds.xCoordinates.values
        y_coords = xrds.y.values if 'y' in xrds else xrds.yCoordinates.values
        
        # Get bounding box (west, south, east, north)
        west = float(np.min(x_coords))
        east = float(np.max(x_coords))
        south = float(np.min(y_coords))
        north = float(np.max(y_coords))
        
        # Create a box geometry (WGS-84 projection assumed for the bounding box)
        bbox = box(west, south, east, north)
        
        return bbox
    except Exception as e:
        print(f"Error extracting bounding box from {h5file}: {str(e)}")
        return None


def get_frameid_from_geojson(bbox, band, trk, path, frame, polar, output_dir=None):
    """
    Create a GeoJSON file containing the frame boundary
    
    Args:
        bbox: shapely Polygon with frame boundaries
        band: Band identifier ('L', 'S', or 'G')
        trk: Track number
        path: Path number
        frame: Frame number
        polar: Polarization ('A' or 'D' or 'HH/VV')
        output_dir: Directory to save the GeoJSON file (default uses metadata dir)
    
    Returns:
        Path to the saved GeoJSON file
    """
    # Determine output directory
    if output_dir is None:
        output_dir = os.path.dirname(os.path.dirname(h5file))
    
    # Create metadata directory if needed
    metadata_dir = os.path.join(output_dir, 'metadata')
    if not os.path.exists(metadata_dir):
        os.makedirs(metadata_dir)
    
    # Create frame ID string for filename
    frame_id = f"NISAR.{trk}.{frame}.{polar}"
    
    # Create GeoDataFrame
    gdf = gpd.GeoDataFrame(
        {'frame': [frame_id], 'band': [band], 'path': [path], 'frame_num': [frame]},
        geometry=[bbox],
        crs='EPSG:4326'
    )
    
    # Save as GeoJSON
    geojson_path = os.path.join(metadata_dir, f"{frame_id}.geojson")
    gdf.to_file(geojson_path, driver='GeoJSON')
    
    # Also save as WKT text file
    wkt_path = os.path.join(metadata_dir, f"{frame_id}.wkt")
    with open(wkt_path, 'w') as f:
        f.write(str(bbox))
    
    print(f"Saved frame boundary to: {geojson_path}")
    print(f"Saved frame WKT to: {wkt_path}")
    
    return geojson_path


def create_nisar_frame_polygons(h5file_list, output_dir=None, output_format='geojson'):
    """
    Process a list of NISAR files and create polygon boundary files for each frame
    
    Args:
        h5file_list: List of NISAR GUNW or GSLC file paths
        output_dir: Output directory (default uses metadata directory of each file)
        output_format: Output format - 'geojson' or 'wkt' or both
    
    Returns:
        Dictionary of frame_id -> polygon geometry for frames that were successfully processed
    """
    results = {}
    
    for h5file in h5file_list:
        print(f"\nProcessing: {h5file}")
        
        # Extract frame info from filename
        try:
            band, trk, path, frame, polar = get_nisar_frame_from_file(h5file)
            frame_id = f"NISAR.{trk}.{frame}.{polar}"
            
            # Extract bounding box
            bbox = extract_frame_bbox_from_gunw(h5file)
            
            if bbox is None:
                print(f"  Warning: Could not extract bounding box for {frame_id}")
                continue
            
            # Save polygon files
            if output_dir:
                output_subdir = output_dir
            else:
                # Use metadata directory of the main file directory
                base_dir = os.path.dirname(os.path.dirname(h5file))
                output_subdir = os.path.join(base_dir, 'metadata')
            
            # Create output GeoDataFrame
            gdf = gpd.GeoDataFrame(
                {'frame': [frame_id], 'band': [band], 'path': [path], 'frame_num': [frame]},
                geometry=[bbox],
                crs='EPSG:4326'
            )
            
            # Save as requested format(s)
            metadata_dir = os.path.join(output_subdir, 'metadata')
            if not os.path.exists(metadata_dir):
                os.makedirs(metadata_dir)
            
            frame_id_clean = f"NISAR_{trk}_{frame}_{polar}"
            geojson_path = os.path.join(metadata_dir, f"{frame_id_clean}.geojson")
            wkt_path = os.path.join(metadata_dir, f"{frame_id_clean}.wkt")
            
            if output_format == 'geojson' or output_format == 'both':
                gdf.to_file(geojson_path, driver='GeoJSON')
                print(f"  Saved GeoJSON to: {geojson_path}")
            
            if output_format == 'wkt' or output_format == 'both':
                with open(wkt_path, 'w') as f:
                    f.write(str(bbox))
                print(f"  Saved WKT to: {wkt_path}")
            
            results[frame_id_clean] = bbox
            
        except Exception as e:
            print(f"  Error processing {os.path.basename(h5file)}: {str(e)}")
            continue
    
    print(f"\nSuccessfully processed {len(results)} frame(s)")
    return results


def export_nisar_metadata_with_polygons(h5file, metafile):
    """
    Extension of export_licsar_metadata to also create polygon boundary files
    
    This function should be called from export_licsar_metadata after metadata.txt is created
    
    Args:
        h5file: Path to NISAR file (GUNW or GSLC)
        metafile: Path to metadata.txt file (will be in same directory)
    """
    # Get the file directory
    file_dir = os.path.dirname(h5file)
    meta_dir = os.path.dirname(metafile)
    
    # Extract frame info
    band, trk, path, frame, polar = get_nisar_frame_from_file(h5file)
    frame_id = f"NISAR_{trk}_{frame}_{polar}"
    
    # Extract bounding box
    bbox = extract_frame_bbox_from_gunw(h5file)
    
    if bbox is None:
        print(f"Warning: Could not extract polygon for {frame_id}")
        return False
    
    # Save as WKT in the metadata directory
    wkt_path = os.path.join(meta_dir, f"{frame_id}.wkt")
    with open(wkt_path, 'w') as f:
        f.write(str(bbox))
    
    print(f"Added polygon boundary to metadata for {frame_id}")
    return True


# For testing
if __name__ == "__main__":
    # Example usage
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python nisardata_polygon.py <h5file> [output_dir]")
        sys.exit(1)
    
    h5file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    # Extract and save polygon
    bbox = extract_frame_bbox_from_gunw(h5file)
    
    if bbox:
        frame_id = get_nisar_frame_from_file(h5file)[0:3]  # band, trk, frame
        frame_id = f"NISAR_{frame_id[1]}{frame_id[2]}_{frame_id[3]}_{frame_id[4]}"
        
        create_nisar_frame_polygons([h5file], output_dir, 'both')
    else:
        print("Failed to extract bounding box")