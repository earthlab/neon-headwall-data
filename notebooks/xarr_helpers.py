import rasterio as rio
import os,sys
from glob import glob
import xarray as xr
from matplotlib import pyplot as plt
import numpy as np
import dask as da
import h5py as h5
import geopandas as gpd
import pandas as pd
from spectral.io import envi
import matplotlib
import dask.array as da
from rasterio import features
from affine import Affine
import warnings
from shapely.geometry.point import Point

# enable fiona KML
import fiona
fiona.drvsupport.supported_drivers['kml'] = 'rw' # enable KML support which is disabled by default
fiona.drvsupport.supported_drivers['KML'] = 'rw' # enable KML support which is disabled by default

# basemap
import contextily as ctx

def add_basemap(ax, zoom, url='http://tile.stamen.com/terrain/tileZ/tileX/tileY.png'):
    xmin, xmax, ymin, ymax = ax.axis()
    basemap, extent = ctx.bounds2img(xmin, ymin, xmax, ymax, zoom=zoom, url=url)
    ax.imshow(basemap, extent=extent, interpolation='bilinear')
    # restore original x/y limits
    ax.axis((xmin, xmax, ymin, ymax))
    
## define some functions for creating xarray datasets out of the files
def NEON_create_rad_xarr_from_h5_file(h5file, nid= 'R10C', nodata=-9999):
    
    # Read H5 file
    f = h5.File(h5file, "r")
    
    # spectral
    wavelength = f[nid]['Radiance']['Metadata']['Spectral_Data']['Wavelength'][:]
    fwhm = f[nid]['Radiance']['Metadata']['Spectral_Data']['FWHM'][:]
    
    # CRS
    crs_str = f[nid]['Radiance']['Metadata']['Coordinate_System']['Coordinate_System_String'].value
    crs_epsg = f[nid]['Radiance']['Metadata']['Coordinate_System']['EPSG Code'].value
    crs_mapinfo = f[nid]['Radiance']['Metadata']['Coordinate_System']['Map_Info'].value
    crs_proj4 = f[nid]['Radiance']['Metadata']['Coordinate_System']['Proj4'].value
    
    #arr = f[nid]['Radiance']['Radiance_Data'][:]
    arr = da.from_array(f[nid]['Radiance']['Radiance_Data'], chunks=(256,256,256))
    
    mapinfo_list = [a.strip() for a in str(crs_mapinfo).split(',')]
    mapinfo = [float(a) for a in mapinfo_list[1:7]]
    mapinfo
    pix_size = mapinfo[0]
    x = np.arange(mapinfo[2], mapinfo[2] + pix_size*arr.shape[1], pix_size)
    y = np.arange(mapinfo[3], mapinfo[3] - pix_size* arr.shape[0], -pix_size)
    
    xr_cube = xr.DataArray(arr, {'y': y, 'x': x, 'bands': wavelength}, dims=['y', 'x', 'bands'])
    xr_cube_ma = xr_cube.where(xr_cube != -9999)
    
    return x, y, xr_cube_ma

def NEON_create_refl_xarr_from_h5_file(h5file, nid='R10C', nodata=-9999):
    
    # Read H5 file
    f = h5.File(h5file, "r")
    
    # spectral
    wavelength = f[nid]['Reflectance']['Metadata']['Spectral_Data']['Wavelength'][:]
    fwhm = f[nid]['Reflectance']['Metadata']['Spectral_Data']['FWHM'][:]
    
    # CRS
    crs_str = f[nid]['Reflectance']['Metadata']['Coordinate_System']['Coordinate_System_String'].value
    crs_epsg = f[nid]['Reflectance']['Metadata']['Coordinate_System']['EPSG Code'].value
    crs_mapinfo = f[nid]['Reflectance']['Metadata']['Coordinate_System']['Map_Info'].value
    crs_proj4 = f[nid]['Reflectance']['Metadata']['Coordinate_System']['Proj4'].value
    
    #arr = f[nid]['Reflectance']['Reflectance_Data'][:]
    arr = da.from_array(f[nid]['Reflectance']['Reflectance_Data'], chunks=(256,256,256))
    sf = f[nid]['Reflectance']['Reflectance_Data'].attrs['Scale_Factor']
    
    mapinfo_list = [a.strip() for a in str(crs_mapinfo).split(',')]
    mapinfo = [float(a) for a in mapinfo_list[1:7]]
    mapinfo
    pix_size = mapinfo[0]
    x = np.arange(mapinfo[2], mapinfo[2] + pix_size*arr.shape[1], pix_size)
    y = np.arange(mapinfo[3], mapinfo[3] - pix_size* arr.shape[0], -pix_size)
    
    xr_cube = xr.DataArray(arr, {'y': y, 'x': x, 'bands': wavelength}, dims=['y', 'x', 'bands'])
    xr_cube_ma = xr_cube.where(xr_cube != -9999)
    
    return x, y, xr_cube_ma/sf
    
def transform_from_latlon(lat, lon):
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    trans = Affine.translation(lon[0], lat[0])
    scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0])
    return trans * scale

def rasterize(shapes, coords, fill=np.nan, **kwargs):
    """Rasterize a list of (geometry, fill_value) tuples onto the given
    xray coordinates. This only works for 1d latitude and longitude
    arrays.
    """
    transform = transform_from_latlon(coords['y'], coords['x'])
    out_shape = (len(coords['y']), len(coords['x']))
    raster = features.rasterize(shapes, out_shape=out_shape,
                                fill=fill, transform=transform,
                                dtype=float, **kwargs)
    return xr.DataArray(raster, coords=coords, dims=('y', 'x'))

def extract_from_headwall(vnir, swir, geodf=None):
    '''Params:
    vnir: list or tuple
        must contain (xarray, x coordinates, y coordinates)
    
    swir: list or tuple
        must contain (xarray, x coordinates, y coordinates)
        
    geom: Shapely geometry
        geometry for which to perform the extraction
        
    Returns: tuple of numpy arrays
        contains column vectors of the spectra for the provided geometry: VNIR, SWIR, and combined    
    
    '''
    
    # check to make sure a geometry is provided
    if geodf is None:
        raise ValueError('Please provide a valid geometry')
        
    
    # parse the inputs
    xarr_vn, x_vnir, y_vnir = vnir
    xarr_sw, x_swir, y_swir = swir
    
    # create mask datasets for the VNIR and SWIR data separately
    ds_vnir = xr.Dataset(coords={'y':y_vnir, 'x':x_vnir})
    shapes = [(shape, n) for n, shape in enumerate(geodf.geometry)]
    ds_vnir['aoi'] = rasterize(shapes, ds_vnir.coords)
    ds_vnir['aoi'] = ds_vnir.aoi + 1

    ds_swir = xr.Dataset(coords={'y':y_swir, 'x':x_swir})
    shapes = [(shape, n) for n, shape in enumerate(geodf.geometry)]
    ds_swir['aoi'] = rasterize(shapes, ds_swir.coords)
    ds_swir['aoi'] = ds_swir.aoi + 1
    
    # apply the mask to the data cube
    example_vnir = ds_vnir.aoi * xarr_vn 
    example_swir = ds_swir.aoi * xarr_sw
    
    # get the valid y and x coordinates, then reduce to unique values
    val_y, val_x = np.where(ds_vnir.aoi==1)
    u_y = np.unique(val_y)
    u_x = np.unique(val_x)
    ex_vnir = example_vnir.sel(y=y_vnir[u_y], x=x_vnir[u_x])
   
    
    val_y, val_x = np.where(ds_swir.aoi==1)
    u_y = np.unique(val_y)
    u_x = np.unique(val_x)
    ex_swir = example_swir.sel(y=y_swir[u_y], x=x_swir[u_x])
    
    print(ex_vnir.shape, ex_swir.shape)
    # shapes may be different....
    if ex_vnir.shape != ex_swir.shape:
        warnings.warn('extracted shapes are not equal, shaving some off...')
        
        min_shape = (min(ex_vnir.shape[0], ex_swir.shape[0]), min(ex_vnir.shape[1], ex_swir.shape[1]), min(ex_vnir.shape[2], ex_swir.shape[2]))
        ex_swir = ex_swir[:min_shape[0], :min_shape[1], :min_shape[2]]
        ex_vnir = ex_vnir[:min_shape[0], :min_shape[1], :min_shape[2]]
        
    print(ex_vnir.shape, ex_swir.shape)
    
    # concatenate the data
    full_ex = np.vstack((ex_vnir.values.reshape(-1, ex_vnir.shape[-1]).T, ex_swir.values.reshape(-1, ex_swir.shape[-1]).T))
    
    
    # concatenate the wavelength vectors
    full_wav = np.concatenate((ex_vnir.coords['wavelength'].values, ex_swir.coords['wavelength'].values))
    full_wav = np.unique(full_wav)
    
    return ex_vnir, ex_swir, (full_wav, full_ex)

def extract_from_headwall_ENVI(vnir, swir, geodf=None, w_cutoff=0):
    '''Params:
    vnir: list or tuple
        must contain (xarray, x coordinates, y coordinates)
    
    swir: list or tuple
        must contain (xarray, x coordinates, y coordinates)
        
    geom: Shapely geometry
        geometry for which to perform the extraction
        
    w_cutoff: float
        wavelength at which to cut off the VNIR spectra
        
    Returns: tuple of numpy arrays
        contains column vectors of the spectra for the provided geometry: VNIR, SWIR, and combined    
    
    '''
    
    # check to make sure a geometry is provided
    if geodf is None:
        raise ValueError('Please provide a valid geometry')
        
    
    # parse the inputs
    xarr_vn, x_vnir, y_vnir = vnir
    xarr_sw, x_swir, y_swir = swir
    
    # create mask datasets for the VNIR and SWIR data separately
    ds_vnir = xr.Dataset(coords={'y':y_vnir, 'x':x_vnir})
    shapes = [(shape, n) for n, shape in enumerate(geodf.geometry)]
    ds_vnir['aoi'] = rasterize(shapes, ds_vnir.coords)
    ds_vnir['aoi'] = ds_vnir.aoi + 1

    ds_swir = xr.Dataset(coords={'y':y_swir, 'x':x_swir})
    shapes = [(shape, n) for n, shape in enumerate(geodf.geometry)]
    ds_swir['aoi'] = rasterize(shapes, ds_swir.coords)
    ds_swir['aoi'] = ds_swir.aoi + 1
    
    # apply the mask to the data cube
    example_vnir = ds_vnir.aoi * xarr_vn 
    example_swir = ds_swir.aoi * xarr_sw
    
    # get the valid y and x coordinates, then reduce to unique values
    val_y, val_x = np.where(ds_vnir.aoi==1)
    u_y = np.unique(val_y)
    u_x = np.unique(val_x)
    ex_vnir = example_vnir.sel(y=y_vnir[u_y], x=x_vnir[u_x])
   
    
    val_y, val_x = np.where(ds_swir.aoi==1)
    u_y = np.unique(val_y)
    u_x = np.unique(val_x)
    ex_swir = example_swir.sel(y=y_swir[u_y], x=x_swir[u_x])
    
    print(ex_vnir.shape, ex_swir.shape)
    # shapes may be different....
    if ex_vnir.shape != ex_swir.shape:
        warnings.warn('extracted shapes are not equal, shaving some off...')
        
        min_shape = (min(ex_vnir.shape[0], ex_swir.shape[0]), min(ex_vnir.shape[1], ex_swir.shape[1]), min(ex_vnir.shape[2], ex_swir.shape[2]))
        ex_swir = ex_swir[:min_shape[0], :min_shape[1], :min_shape[2]]
        ex_vnir = ex_vnir[:min_shape[0], :min_shape[1], :min_shape[2]]
        
        print(ex_vnir.shape, ex_swir.shape, min_shape)
        
    # cut off the VNIR if provided
    if w_cutoff > 0:
        b_cutoff = np.where(ex_vnir.coords['wavelength'] <= w_cutoff)[0][-1] + 1 # +1 due to 1 based indexing on band
        ex_vnir = ex_vnir.sel(band=slice(0, b_cutoff))
    
    
    # concatenate the data
    full_ex = np.vstack((ex_vnir.values.reshape(-1, ex_vnir.shape[-1]).T, ex_swir.values.reshape(-1, ex_swir.shape[-1]).T))
    
    
    # concatenate the wavelength vectors
    full_wav = np.concatenate((ex_vnir.coords['wavelength'].values, ex_swir.coords['wavelength'].values))
    full_wav = np.unique(full_wav)
    
    return ex_vnir, ex_swir, (full_wav, full_ex)

def extract_from_NEON(hsi, geodf=None):
    '''Params:
    hsi: list or tuple
        must contain (xarray, x coordinates, y coordinates)
        
    geom: Shapely geometry
        geometry for which to perform the extraction
        
    Returns: tuple of numpy arrays
        contains column vectors of the spectra for the provided geometry: VNIR, SWIR, and combined    
    
    '''
    
    # check to make sure a geometry is provided
    if geodf is None:
        raise ValueError('Please provide a valid geometry')
        
        
    # parse the inputs
    xarr_n, x_n, y_n = hsi
    
    # NEON's data is in UTM, hence to to_crs() call
    ds_neon = xr.Dataset(coords={'y':y_n, 'x':x_n})
    shapes = [(shape, n) for n, shape in enumerate(geodf.geometry)]
    ds_neon['aoi'] = rasterize(shapes, ds_neon.coords)
    ds_neon['aoi'] = ds_neon.aoi + 1

    example_neon = ds_neon.aoi * xarr_n
    
    val_y, val_x = np.where(ds_neon.aoi>=1)
    u_y = np.unique(val_y)
    u_x = np.unique(val_x)
    ex_neon = example_neon.sel(y=y_n[u_y], x=x_n[u_x])

    full_neon = ex_neon.values.reshape(-1, ex_neon.shape[-1]).T
    neon_wav = ex_neon.coords['bands'].values
    
    return ex_neon, (neon_wav, full_neon)

def extract_from_NEON_ENVI(hsi, geodf=None):
    '''Params:
    hsi: list or tuple
        must contain (xarray, x coordinates, y coordinates)
        
    geom: Shapely geometry
        geometry for which to perform the extraction
        
    Returns: tuple of numpy arrays
        contains column vectors of the spectra for the provided geometry: VNIR, SWIR, and combined    
    
    '''
    
    # check to make sure a geometry is provided
    if geodf is None:
        raise ValueError('Please provide a valid geometry')
        
        
    # parse the inputs
    xarr_n, x_n, y_n = hsi
    
    # NEON's data is in UTM, hence to to_crs() call
    ds_neon = xr.Dataset(coords={'y':y_n, 'x':x_n})
    shapes = [(shape, n) for n, shape in enumerate(geodf.geometry)]
    ds_neon['aoi'] = rasterize(shapes, ds_neon.coords)
    ds_neon['aoi'] = ds_neon.aoi + 1

    example_neon = ds_neon.aoi * xarr_n
    
    val_y, val_x = np.where(ds_neon.aoi==1)
    u_y = np.unique(val_y)
    u_x = np.unique(val_x)
    ex_neon = example_neon.sel(y=y_n[u_y], x=x_n[u_x])

    full_neon = ex_neon.values.reshape(-1, ex_neon.shape[-1]).T
    neon_wav = ex_neon.coords['wavelength'].values
    
    return ex_neon, (neon_wav, full_neon)
        
def bin_data(hw_spectra, hw_wavelengths, neon_spectra, neon_wavelengths, res=200):

    # assign inputs
    full_wav = hw_wavelengths
    full_ex = hw_spectra
    
    neon_wav = neon_wavelengths
    full_neon_ls = neon_spectra
    
    
    #res = 200
    neon_std_ls, hw_std_ls = [],[]
    neon_mean_ls, hw_mean_ls = [],[]
    band_ranges = []
    for bc in np.arange(full_wav[0]+res/2, full_wav[-1], res):
        w_min = bc - res/2
        w_max = bc + res/2
        band_ranges.append((w_min, bc, w_max))
        #print(w_min,bc, w_max)

        # get min bands for headwall
        h_b_min = np.where(full_wav >= w_min)[0][0] + 1 # +1 due to 1 based indexing on band
        h_b_max = np.where(full_wav <= w_max)[0][-1] + 1 # +1 due to 1 based indexing on band

        # record std and mean
        hw_std = np.nanstd(full_ex[h_b_min:h_b_max,:])
        hw_std_ls.append(hw_std)
        hw_mean = np.nanmean(full_ex[h_b_min:h_b_max,:])
        hw_mean_ls.append(hw_mean)

        # get min bands for NEON
        n_b_min = np.where(neon_wav >= w_min)[0][0] + 1 # +1 due to 1 based indexing on band
        n_b_max = np.where(neon_wav <= w_max)[0][-1] + 1 # +1 due to 1 based indexing on band

        # record std for each NEON group
        temp_, temp__ = [],[]
        for full_neon in full_neon_ls:
            neon_std = np.nanstd(full_neon[n_b_min:n_b_max,:])
            temp_.append(neon_std)

            neon_mean = np.nanmean(full_neon[n_b_min:n_b_max,:])
            temp__.append(neon_mean)

        neon_std_ls.append(temp_)
        neon_mean_ls.append(temp__)
        
    arr_std = np.hstack((np.expand_dims(np.array(hw_std_ls),axis=1), np.array(neon_std_ls)))
    arr_mean = np.hstack((np.expand_dims(np.array(hw_mean_ls),axis=1), np.array(neon_mean_ls)))
    
    return arr_mean, arr_std

def bin_data_groups(hw_spectra, hw_wavelengths, neon_spectra, neon_wavelengths, res=200):

    # assign inputs
    full_wav = hw_wavelengths
    full_ex_ls = hw_spectra   # list
    
    neon_wav = neon_wavelengths
    full_neon_ls = neon_spectra  # list
    
    
    #res = 200
    neon_std_ls, hw_std_ls = [],[]
    neon_mean_ls, hw_mean_ls = [],[]
    band_ranges = []
    for bc in np.arange(full_wav[0]+res/2, full_wav[-1], res):
        w_min = bc - res/2
        w_max = bc + res/2
        band_ranges.append((w_min, bc, w_max))
        #print(w_min,bc, w_max)

        # get min bands for headwall
        h_b_min = np.where(full_wav >= w_min)[0][0] + 1 # +1 due to 1 based indexing on band
        h_b_max = np.where(full_wav <= w_max)[0][-1] + 1 # +1 due to 1 based indexing on band

        # record std and mean for each HW group
        hw_temp_, hw_temp__ = [],[]
        for full_ex in full_ex_ls:
            
            if len(full_ex.shape) < 2:
                full_ex = full_ex[:, np.newaxis]
                
            hw_std = np.nanstd(full_ex[h_b_min:h_b_max,:])
            hw_temp_.append(hw_std)
            
            hw_mean = np.nanmean(full_ex[h_b_min:h_b_max,:])
            hw_temp__.append(hw_mean)
            
        hw_std_ls.append(hw_temp_)            
        hw_mean_ls.append(hw_temp__)

        # get min bands for NEON
        n_b_min = np.where(neon_wav >= w_min)[0][0] + 1 # +1 due to 1 based indexing on band
        n_b_max = np.where(neon_wav <= w_max)[0][-1] + 1 # +1 due to 1 based indexing on band

        # record std for each NEON group
        temp_, temp__ = [],[]
        for full_neon in full_neon_ls:
            
            if len(full_neon.shape) < 2:
                full_neon = full_neon[:, np.newaxis]
                
            neon_std = np.nanstd(full_neon[n_b_min:n_b_max,:])
            temp_.append(neon_std)

            neon_mean = np.nanmean(full_neon[n_b_min:n_b_max,:])
            temp__.append(neon_mean)

        neon_std_ls.append(temp_)
        neon_mean_ls.append(temp__)
        
    arr_std = np.hstack((np.expand_dims(np.array(hw_std_ls),axis=1).squeeze(), np.array(neon_std_ls)))
    arr_mean = np.hstack((np.expand_dims(np.array(hw_mean_ls),axis=1).squeeze(), np.array(neon_mean_ls)))
    
    return arr_mean, arr_std
