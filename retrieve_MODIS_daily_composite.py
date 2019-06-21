import os, sys
# import urllib2 # use this for python 2x
import urllib.request as urllib2
import subprocess
from osgeo import gdal
import datetime
# import HTMLParser
import re
import lxml.html
from io import BytesIO
from lxml import etree
import numpy as np
import math
import geopandas as gpd
from pprint import pprint
import rasterio as rio

def retrieveMODISlistings(url):

    tree = etree.HTML(urllib2.urlopen(url_base).read())

    # data = urllib2.urlopen(url_base).read()
    # dom = lxml.html.parse(BytesIO(data))
    # xpatheval=etree.XPathDocumentEvaluator(dom)
    # xpatheval('//*[@id="ftp-directory-list"]/tbody/tr[2]/td[1]/a')
    
    # the xpath argument retrieved from F12 mode in chrome
    table = tree.xpath('//*[@id="ftp-directory-list"]')[0]
    
    # grab the hdf files (tested and probably only works for the MODIS site.)
    hdf_files = [row[0][0].items()[0][1] for row in table[2:]]
        
    return hdf_files
    
def sortByTOD(flist, timeofday, num=1):

    hourmin = int(timeofday[:4])    
    times = [int(f.split('/')[-1].split('.')[2]) for f in flist]
    
    times_arr = np.abs(np.array(times) - hourmin)
    inds = list(np.argsort(times_arr))
    
    sorted = [flist[i] for i in inds]
    
    return sorted[:num]
    
    

def modisInfoFromAster(file_name):

    aster = gdal.Open(file_name)
    aster_sds = aster.GetSubDatasets()
    meta = aster.GetMetadata()
    # aster.Close()
    # aster = None

    ## extract the fields (they will be strings)
    tod = meta['TIMEOFDAY']
    calendar_date = meta['CALENDARDATE']
    vnir_view_Z = meta['POINTINGANGLE.1']
    swir_view_Z = meta['POINTINGANGLE.2']
    solar_d = meta['SOLARDIRECTION'] # 2 elements... azimuth and zenith (i think)
    
    return calendar_date, tod
    
def viewInfoFromAster(file_name):

    aster = gdal.Open(file_name)
    meta = aster.GetMetadata()
    # aster.Close()
    # aster = None

    ## extract the fields (they will be strings)
    tod = meta['TIMEOFDAY']
    calendar_date = meta['CALENDARDATE']
    vnir_view_Z = meta['POINTINGANGLE.1']
    swir_view_Z = meta['POINTINGANGLE.2']
    solar_d = meta['SOLARDIRECTION'] # 2 elements... azimuth and zenith (i think)
    
    solar_d = [float(s) for s in solar_d.split(',')]
    vnir_view_Z = float(vnir_view_Z)
    swir_view_Z = float(swir_view_Z)    
    
    return vnir_view_Z, swir_view_Z, solar_d

def bboxInfoFromAster(file_name):

    aster = gdal.Open(file_name)
    meta = aster.GetMetadata()
    # aster.Close()
    # aster = None

    ## extract the fields (they will be strings)
    n = float(meta['NORTHBOUNDINGCOORDINATE'])
    s = float(meta['SOUTHBOUNDINGCOORDINATE'])
    e = float(meta['EASTBOUNDINGCOORDINATE'])
    w = float(meta['WESTBOUNDINGCOORDINATE'])    
    
    return n,s,e,w

def summarizeMODIS_mod04(hdf, aster_bbox):

    mod = gdal.Open(hdf)
    mod_sds = mod.GetSubDatasets()
    
    for sd in mod_sds:
    
        if 'Longitude' in sd[1]:
            longitude_sd = sd[0]
        
        if 'Latitude' in sd[1]:
            latitude_sd = sd[0]
        
        if 'Deep_Blue_Aerosol_Optical_Depth_550_Land' == sd[0].split(':')[-1]:
            aod_sd = sd[0]
            
    # 'close' the file
    mod = None
    
    ## should now have the file paths to open these arrays and sample within the aster_bbox
    lon_arr = gdal.Open(longitude_sd).ReadAsArray()
    lat_arr = gdal.Open(latitude_sd).ReadAsArray()
    aod_arr = gdal.Open(aod_sd).ReadAsArray()
    
    # get the scale and offset
    scale_factor_aod = float(gdal.Open(aod_sd).GetMetadata()['scale_factor'])
    add_offset_aod = float(gdal.Open(aod_sd).GetMetadata()['add_offset'])
    
    print(scale_factor_aod, add_offset_aod)
    
    # create the logical array for the AOI
    lat_max, lat_min, lon_min, lon_max = aster_bbox
    sample_arr = np.where((lon_arr > lon_min) & (lon_arr < lon_max) & (lat_arr > lat_min) & (lat_arr < lat_max))[0]
    
    # sample the data
    temp_arr = aod_arr[sample_arr]
    temp_arr = (temp_arr[temp_arr > -9999] - add_offset_aod) * scale_factor_aod
    aod_mean = np.mean(temp_arr)
            
    return temp_arr, aod_mean
    
    
def summarizeMODIS_mod05(hdf, aster_bbox):

    mod = gdal.Open(hdf)
    mod_sds = mod.GetSubDatasets()
    
    for sd in mod_sds:
    
        if 'Longitude' in sd[1]:
            longitude_sd = sd[0]
        
        if 'Latitude' in sd[1]:
            latitude_sd = sd[0]
        
        # use the Water_Vapor_Near_Infrared dataset to avoid problems (known issue)
        if 'Water_Vapor_Near_Infrared' == sd[0].split(':')[-1]:
            wv_sd = sd[0]
            
    # 'close' the file
    mod = None
    
    ## should now have the file paths to open these arrays and sample within the aster_bbox
    lon_arr = gdal.Open(longitude_sd).ReadAsArray()
    lat_arr = gdal.Open(latitude_sd).ReadAsArray()
    wv_arr = gdal.Open(wv_sd).ReadAsArray()
    
    # get the scale and offset
    scale_factor_wv = float(gdal.Open(wv_sd).GetMetadata()['scale_factor'])
    add_offset_wv = float(gdal.Open(wv_sd).GetMetadata()['add_offset'])
    
    print(scale_factor_wv, add_offset_wv)
    
    # create the logical array for the AOI
    lat_max, lat_min, lon_min, lon_max = aster_bbox
    sample_arr = np.where((lon_arr > lon_min) & (lon_arr < lon_max) & (lat_arr > lat_min) & (lat_arr < lat_max))[0]
    
    # sample the data
    temp_arr = wv_arr[sample_arr]
    temp_arr = (temp_arr[temp_arr > -9999] - add_offset_wv) * scale_factor_wv
    wv_mean = np.mean(temp_arr)
            
    return temp_arr, wv_mean
    
    
def summarizeMODIS_mod07(hdf, aster_bbox):

    mod = gdal.Open(hdf)
    mod_sds = mod.GetSubDatasets()
    
    for sd in mod_sds:
    
        if 'Longitude' in sd[1]:
            longitude_sd = sd[0]
        
        if 'Latitude' in sd[1]:
            latitude_sd = sd[0]
        
        if 'Total_Ozone' == sd[0].split(':')[-1]:
            oz_sd = sd[0]
            
        if 'Water_Vapor' == sd[0].split(':')[-1]:
            wv_sd = sd[0]
            
    # 'close' the file
    mod = None
    
    ## should now have the file paths to open these arrays and sample within the aster_bbox
    lon_arr = gdal.Open(longitude_sd).ReadAsArray()
    lat_arr = gdal.Open(latitude_sd).ReadAsArray()
    oz_arr = gdal.Open(oz_sd).ReadAsArray()
    wv_arr = gdal.Open(wv_sd).ReadAsArray()
    
    # get the scale and offset
    scale_factor_oz = float(gdal.Open(oz_sd).GetMetadata()['scale_factor'])
    add_offset_oz = float(gdal.Open(oz_sd).GetMetadata()['add_offset'])
    
    scale_factor_wv = float(gdal.Open(wv_sd).GetMetadata()['scale_factor'])
    add_offset_wv = float(gdal.Open(wv_sd).GetMetadata()['add_offset'])
    
    print(scale_factor_oz, add_offset_oz, scale_factor_wv, add_offset_wv)
    
    # get the logical array for the AOI
    lat_max, lat_min, lon_min, lon_max = aster_bbox
    sample_arr = np.where((lon_arr > lon_min) & (lon_arr < lon_max) & (lat_arr > lat_min) & (lat_arr < lat_max))[0]
    
    # sample the values
    temp_arr_oz = oz_arr[sample_arr]
    temp_arr_oz = (temp_arr_oz[temp_arr_oz > -9999] - add_offset_oz) * scale_factor_oz
    
    temp_arr_wv = wv_arr[sample_arr]
    temp_arr_wv = (temp_arr_wv[temp_arr_wv > -9999] - add_offset_wv) * scale_factor_wv
    
    oz_mean = np.mean(temp_arr_oz)
    wv_mean = np.mean(temp_arr_wv)
            
    return temp_arr_oz, temp_arr_wv, oz_mean, wv_mean
        

def summarizeMODIS_mod08_d3(hdf, aster_bbox, variable_name, pt):

    mod = gdal.Open(hdf)
    mod_sds = mod.GetSubDatasets()
    
    for sd in mod_sds:
    
        # use the Water_Vapor_Near_Infrared dataset to avoid problems (known issue)
        if variable_name == sd[0].split(':')[-1]:
            mod_sd = sd[0]
            
    # 'close' the file
    mod = None
    
    # get the scale and offset for the subdataset
    scale_factor = float(gdal.Open(mod_sd).GetMetadata()['scale_factor'])
    add_offset = float(gdal.Open(mod_sd).GetMetadata()['add_offset'])
    
    # sample the subdataset at XY defined by pt
    samp_vals = []
    with rio.open(mod_sd) as src:
        for val in src.sample([(pt[0], pt[1])]):
            proc_val = (val - add_offset) * scale_factor
            if val == -9999:
                proc_val = [None]
            samp_vals.append(proc_val[0])
    
    return samp_vals
    


# change directory to correct path
dl_dir = r"D:\projects\headwall_neon"
save_dir = r'D:\projects\headwall_neon' # make parameter
os.chdir(dl_dir)

#################################################################################################################
### ASTER METADATA PROCESSING ###
#################################################################################################################

# specify the ASTER file (potentially download some other way... for now, requested through EarthData) 
# make parameter
#aster_file = r'C:\tools\py6s_emulator\6S_emulator\examples\AST_L1T_00309132006180645_20150516042027_58950.hdf'

# get the bounding box
# n,s,e,w = bboxInfoFromAster(aster_file)
# aster_bbox = [n,s,w,e]
bbox = np.array([[-105.24818658828735,40.12792502825585], 
                    [-105.24319767951965,40.12792502825585],
                    [-105.24319767951965,40.13174768054406],
                    [-105.24818658828735,40.13174768054406],
                    [-105.24818658828735,40.12792502825585]])
aster_bbox = [np.max(bbox[:,1]), np.min(bbox[:,1]), np.min(bbox[:,0]), np.max(bbox[:,0])]
aster_bbox = [np.max(bbox[:,1]), np.min(bbox[:,1]), np.max(bbox[:,0]), np.min(bbox[:,0])]

## specify center point, buffer to 4km
from shapely.geometry import Point
from fiona.crs import from_epsg

pt_loc = Point(-105.24561166763304,40.12962309991011)
pt_loc2 = (-105.24561166763304,40.12962309991011)
df = gpd.GeoDataFrame({'geometry':[pt_loc]})
df.crs = from_epsg(4326)
df_proj = df.to_crs(epsg=3857)
bounds = df_proj.buffer(10000).to_crs(epsg=4326).bounds
aster_bbox = [bounds.maxy[0], bounds.miny[0], bounds.minx[0], bounds.maxx[0]]

# get the calendar_Date
#calendar_date, tod = modisInfoFromAster(aster_file)

# convert the calendar_date to year / month / day
# year  = int(calendar_date[0:4])
# month = int(calendar_date[4:6])
# day   = int(calendar_date[6:])
year = 2019
month = 4
day = 9
tod = '1630' # add 6 hours for GMT, was 1030 local time

# convert calendar date to doy
d = datetime.date(year, month, day)
doy = d.timetuple().tm_yday

#################################################################################################################
#################################################################################################################
#################################################################################################################


# info for downloading modis data
# example file https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/
#                       MOD05_L2/2006/066/MOD05_L2.A2006066.0000.061.2017263153530.hdf
aerosol_prod = 'MOD04_L2' #daily, MOD--Terra
total_precip = 'MOD05_L2' #daily
atm_profile  = 'MOD07_L2' #daily

daily_prod = 'MOD08_D3'
daily_prod = 'MYD08_D3'


# get the MODIS HDF files
url_base = 'https://ladsweb.modaps.eosdis.nasa.gov/archive/allData/61/{0}/{1}/{2:0>3}/'.format(daily_prod, year, doy)
potential_files = retrieveMODISlistings(url_base)
hdf_file = os.path.basename(potential_files[0])
url = url_base + hdf_file
full_fname = os.path.join(save_dir, hdf_file)
if not os.path.exists(full_fname):
    filedata = urllib2.urlopen(url)
    datatowrite=filedata.read()
    with open(full_fname, 'wb') as f:
        f.write(datatowrite)

# sample the data
vals = []
prods = ('Water_Vapor_Near_Infrared_Clear_Mean', 
         'Water_Vapor_Near_Infrared_Cloud_Mean',
         'Total_Ozone_Mean',
         'Deep_Blue_Aerosol_Optical_Depth_550_Land_Mean')
         
for prod in prods:

    var = summarizeMODIS_mod08_d3(full_fname, aster_bbox, prod, pt_loc2)
    vals.append(var)
    
for v in vals:
    print(v[0])