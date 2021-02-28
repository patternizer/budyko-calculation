#! /usr/bin python

#------------------------------------------------------------------------------
# PROGRAM: budyko-calculation.py
#------------------------------------------------------------------------------
# Version 0.3
# 28 February, 2021
# Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
# michael DOT a DOT taylor AT uea DOT ac DOT uk 
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# IMPORT PYTHON LIBRARIES
#------------------------------------------------------------------------------
# Numerics and dataframe libraries:
import numpy as np
import pandas as pd
import xarray as xr
from datetime import datetime
import calendar as cal
# Plotting libraries:
import matplotlib
import matplotlib.pyplot as plt; plt.close('all')
import matplotlib.cm as cm
from matplotlib import rcParams
from matplotlib import colors as mcol
from matplotlib.cm import ScalarMappable
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import matplotlib.path as mpath
from matplotlib.collections import PolyCollection
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
import cmocean
# Mapping libraries:
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cf
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
matplotlib.rcParams['text.usetex'] = True
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Avant Garde', 'Lucida Grande', 'Verdana', 'DejaVu Sans' ]
plt.rc('savefig',facecolor='white')
plt.rc('axes',edgecolor='black')
plt.rc('xtick',color='black')
plt.rc('ytick',color='black')
plt.rc('axes',labelcolor='black')
plt.rc('axes',facecolor='white')

#----------------------------------------------------------------------------
# SETTINGS
#----------------------------------------------------------------------------

fontsize = 16
co2_file = 'monthly_in_situ_co2_mlo.csv'
ice_threshold = 0.1
use_polar = True
use_rolling = True
use_era5 = True
plot_polar_maps = True
plot_versus_time = True
if use_rolling == True:
    avestr = 'ma-12'
else:    
    avestr = 'yearly'
if use_era5 == True:
    reanalysisstr = 'ERA5'
else:    
    reanalysisstr = 'JRA55'
thresholdstr = str(ice_threshold)

# Calculate current time for polar plots

now = datetime.now()
currentmn = str(now.month)
if now.day == 1:
    currentdy = str(cal.monthrange(now.year,now.month-1)[1])
    currentmn = str(now.month-1)
else:
    currentdy = str(now.day-1)
if int(currentdy) < 10:
    currentdy = '0' + currentdy    
currentyr = str(now.year)
if int(currentmn) < 10:
    currentmn = '0' + currentmn
currenttime = currentmn + '_' + str(currentdy) + '_' + currentyr
titletime = currentmn + '/' + str(currentdy) + '/' + currentyr

#----------------------------------------------------------------------------
# LOAD: CO2 monthly timeseries
#----------------------------------------------------------------------------

nheader = 6
f = open(co2_file)
lines = f.readlines()
years = []
months = []
co2_fits = []
co2_fits_sa = [] # seasonally_adjusted
for i in range(nheader,len(lines)):
    words = lines[i].split(',')    
    year = words[0]
    month = words[1]
    co2_fit = words[8]
    co2_fit_sa = words[9].rstrip()
    years.append(year)
    months.append(month)
    co2_fits.append(float(co2_fit.strip()))
    co2_fits_sa.append(float(co2_fit_sa.strip()))
f.close()    
years = np.array(years)
months = np.array(months)
days = np.ones(np.size(years))
co2_fits = np.array(co2_fits)
co2_fits_sa = np.array(co2_fits_sa)
co2_dates = [str(years[i].astype(int)) + '-' + str(months[i].astype(int)).zfill(2) + '-' + str(days[i].astype(int)).zfill(2) for i in range(len(years))]

df = pd.DataFrame({'date':co2_dates,'co2_fits_sa_monthly':co2_fits_sa})
df['co2_fits_sa_monthly'].replace(-99.99, np.nan, inplace=True)

#----------------------------------------------------------------------------
# LOAD: sea ice field: sic1x1 (1=ice, 0=none)
#----------------------------------------------------------------------------

if use_era5 == True:    
    ds = xr.open_dataset('era5_sic_1x1.nc')
else:
    ds = xr.open_dataset('jra55_sic_1x1.nc')
lat = ds.lat
lon = ds.lon 
par = ds.sic
time = ds.time
t = [time[i].values.ravel()[0] for i in range(len(time))]
sic_dates = [str(t[i])[0:4] + '-' + str(t[i])[4:6] + '-' + str(t[i])[6:8] for i in range(len(time))]
lats = np.array(par.lat)
lons = np.array(par.lon)
icemin_avelat = []

N = par.shape[0]
for i in range(N):
    var = par[i,:,:]
    var_region = var.where((var>=ice_threshold)&(var['lat']>0))
    result = np.array(var_region)
    
    # Calculate mean latitudes
    
    latmins = []
    for j in range(len(lons)):
        ice = result[:,j] > 0   
        if np.sum(ice) == 0:
            latmin = np.nan
        else:                
            latmin = lats[ice.nonzero()[0][0]]
        latmins.append(latmin)

    # Require at least 30 longitudinal points to estimate mean lat

    if (np.sum(latmins) <= 30):
        avelat = np.nan
    else:        
        avelat = np.nanmedian(latmins)
    icemin_avelat.append(avelat)

    #----------------------------------------------------------------------------
    # PLOT: sea ice cover (polar view)
    #----------------------------------------------------------------------------
    
    # PLOT: sea-ice cover on 1969-01
    #
    # if use_era5 == True:        
    #    iceland_test_year_idx = 228
    # else:
    #    iceland_test_year_idx = 133
    # if i == iceland_test_year_idx:

    if plot_polar_maps == True:

        figstr = 'sea-ice-' + sic_dates[i] + '-' + reanalysisstr + '-' + avestr + '-' + thresholdstr + '.png'             
        titlestr = 'NH sea-ice covererage on ' + sic_dates[i]  + ': ' + r'sic$\ge$' + thresholdstr + r' N$\ge$30 (for median)'            

        fig  = plt.figure(figsize=(15,10))
        cmap = cmocean.cm.ice
        data = var_region
        lat = np.array(data.lat)
        lon = np.array(data.lon)
        z = var_region[:,:]
        x, y = np.meshgrid(lon, lat) # convert vector lat and lon to 2D
        if use_polar == True:
            p = ccrs.NorthPolarStereo(); threshold = 0
        else:
            p = ccrs.PlateCarree(central_longitude=0); threshold = 0
        ax = plt.axes(projection=p)
        ax.add_feature(cf.LAND, facecolor='linen')
        ax.add_feature(cf.COASTLINE, edgecolor="black", linewidth=0.7)
        ax.add_feature(cf.BORDERS, edgecolor="black", linewidth=0.5)
        g = ccrs.Geodetic()           
        trans = ax.projection.transform_points(g, x, y)
        x0 = trans[:,:,0]
        y0 = trans[:,:,1]
        if use_polar == True:
            ax.set_extent([-180,180,90,30], ccrs.PlateCarree())
            parallels = np.arange(30,86,10)
            meridians = np.arange(-180,180,30)
            gl = ax.gridlines(crs=ccrs.PlateCarree(), xlocs=meridians, ylocs=parallels, linestyle="dotted", linewidth=0.7, color='black', alpha=0.5)
            theta = np.linspace(0, 2*np.pi, 120)
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            center, radius = [0.5, 0.5], 0.5
            circle = mpath.Path(verts * radius + center)
            ax.set_boundary(circle, transform=ax.transAxes)  # circular bounding box (not rectangular)
            plt.contourf(x, y, z, transform=ccrs.PlateCarree(), levels=10, cmap=cmap)   
            cb = plt.colorbar(orientation="vertical", shrink=0.5, pad=0.05)    
            cb.ax.set_title('sic', fontsize=fontsize)        
            cb.ax.tick_params(labelsize=fontsize)
            plt.scatter(x=lon, y=latmins, c='red', s=3, alpha=1.0, transform=ccrs.PlateCarree(), label='minimum')   
            plt.plot(lon, [avelat]*360, color='gold', lw=3, alpha=1.0, transform=ccrs.PlateCarree(), label='median')  
            plt.legend(bbox_to_anchor=(1.0, 1), loc='upper left', fontsize=fontsize, fancybox=True)
        else:
            ax.set_extent([-180, 180, 90, 60], crs=p)    
            gl = ax.gridlines(crs=p, draw_labels=True, linewidth=1, color='lightgrey', alpha=1.0, linestyle='-')
            gl.xlabels_top = False
            gl.ylabels_right = False
            gl.xlines = True
            gl.ylines = True
            gl.xlocator = mticker.FixedLoca2tor([-180,-120,-60,0,60,120,180])
            gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
            gl.xformatter = LONGITUDE_FORMATTER
            gl.yformatter = LATITUDE_FORMATTER 
            gl.xlabel_style = {'size': fontsize}
            gl.ylabel_style = {'size': fontsize}        
            plt.scatter(x=x, y=y, c=z[:,:], s=1, alpha=1.0, transform=p)         
            cb = plt.colorbar(orientation="vertical", shrink=0.5, pad=0.05)    
            cb.ax.set_title('sic', fontsize=fontsize)        
            cb.ax.tick_params(labelsize=fontsize)

        fig.suptitle(titlestr, fontsize=24, color='black', fontweight='bold')
        plt.annotate(r'$\textbf{Data}$' + ': ' + reanalysisstr, xy=(150,90),
             xycoords='figure pixels', color='black', fontsize=fontsize) 
        plt.annotate(r'$\textbf{Graphic}$' + ': Michael Taylor, CRU/UEA (@MichaelTaylorEO)' + ' -- ' + titletime, xy=(150,50),
             xycoords='figure pixels', color='black', fontsize=fontsize,
             bbox=dict(boxstyle="square, pad=0.3", fc='white', edgecolor='black', linewidth=0.2))     
        fig.subplots_adjust(top=0.92)
        plt.savefig(figstr)        
        plt.close(fig)

dg = pd.DataFrame({'date':sic_dates,'icemin_avelat_monthly':icemin_avelat})

# TRIM:  dataframes to common time range

tmin = df.min()[0]
tmax = dg.max()[0]
df_trim = df[(df['date']>=tmin)&(df['date']<=tmax)].reset_index(drop=True)
dg_trim = dg[(dg['date']>=tmin)&(dg['date']<=tmax)].reset_index(drop=True)
yyyy = []
mm = []
dd = []
for i in range(len(df_trim)):    
    y = int(df_trim['date'][i][0:4])
    m = int(df_trim['date'][i][5:7])
    d = int(df_trim['date'][i][8:10])
    yyyy.append(y)
    mm.append(m)
    dd.append(d)   
times = [datetime(yyyy[i],mm[i],dd[i]) for i in range(len(df_trim))]

# CALCULATE: yearly average

if use_rolling == True:
 
    df_trim['co2_fits_sa_yearly'] = df_trim['co2_fits_sa_monthly'].dropna().rolling(12).mean()
    dg_trim['icemin_avelat_yearly'] = dg_trim['icemin_avelat_monthly'].dropna().rolling(12).mean()
    
else:

    df_yearly = []
    for i in range(len(np.unique(years))):
        yearly_ave = [np.nanmean(df_trim[years==np.unique(years)[i]]['co2_fits_sa_monthly'])]*12
        df_yearly += yearly_ave
    df_yearly = df_yearly[0:len(df_trim)]
    df_trim['co2_fits_sa_yearly'] = df_yearly
    
    dg_yearly = []
    for i in range(len(np.unique(years))):
        yearly_ave = [np.nanmean(dg_trim[years==np.unique(years)[i]]['icemin_avelat_monthly'])]*12
        dg_yearly += yearly_ave
    dg_yearly = dg_yearly[0:len(dg_trim)]
    dg_trim['icemin_avelat_yearly'] = dg_yearly

file_str = reanalysisstr + '-' + avestr + '-' + thresholdstr + '.csv'
df_str = 'df'+ '-' + file_str
dg_str = 'dg'+ '-' + file_str
df_trim.to_csv(df_str)
dg_trim.to_csv(dg_str)

#----------------------------------------------------------------------------
# PLOT: sea ice minimum Arctic latitude median vs CO₂ (and temporally)
#----------------------------------------------------------------------------

if plot_versus_time == True:

    figstr = 'budyko-calculation-' + reanalysisstr + '-' + avestr + '-' + thresholdstr + '.png'
    titlestr = r'Sea-ice boundary median latitude versus CO$_{2}$: ' + reanalysisstr + ' (' + avestr + '): ' + r'sic$\ge$' + thresholdstr + r' N$\ge$30 (for median)'
                 
    fig, ax = plt.subplots(figsize=(15,10))          
    plt.step(df_trim['co2_fits_sa_yearly'],dg_trim['icemin_avelat_yearly'], color='black', alpha=1.0)
    plt.ylim(60,90)
    ax.xaxis.grid(True, which='major')      
    ax.yaxis.grid(True, which='major')  
    plt.tick_params(labelsize=fontsize)    
    plt.xlabel(r'CO$_{2}$ [ppmv]', fontsize=fontsize)
    plt.ylabel(r'Median latitude (sea-ice min) [$^{\circ}$N]', fontsize=fontsize)
    plt.title(titlestr, fontsize=fontsize)
    plt.savefig(figstr)
    plt.close(fig)

    figstr = 'budyko-calculation-' + reanalysisstr + '-' + 'co2' + '-' + avestr + '-' + thresholdstr + '.png'
    titlestr = r'CO$_{2}$ with time: ' + avestr

    fig, ax = plt.subplots(figsize=(15,10))          
    plt.step(times,df_trim['co2_fits_sa_yearly'], color='black', alpha=1.0)
    plt.xlim('1955-01-01','2025-01-01')
    ax.xaxis.grid(True, which='major')      
    ax.yaxis.grid(True, which='major')  
    plt.tick_params(labelsize=fontsize)    
    plt.ylabel(r'CO$_{2}$ [ppmv]', fontsize=fontsize)
    plt.xlabel('Time', fontsize=fontsize)
    plt.title(titlestr, fontsize=fontsize)
    plt.savefig(figstr)
    plt.close(fig)

    figstr = 'budyko-calculation-' + reanalysisstr + '-' + 'time' + '-' + avestr + '-' + thresholdstr + '.png'
    titlestr = 'Sea-ice boundary median latitude with time: ' + reanalysisstr + ' (' + avestr + '): ' + r'sic$\ge$' + thresholdstr + r' N$\ge$30 (for median)'

    fig, ax = plt.subplots(figsize=(15,10))          
    plt.step(times,dg_trim['icemin_avelat_yearly'], color='black', alpha=1.0)
    plt.xlim('1955-01-01','2025-01-01')
    plt.ylim(60,90)
    ax.xaxis.grid(True, which='major')      
    ax.yaxis.grid(True, which='major')  
    plt.tick_params(labelsize=fontsize)    
    plt.ylabel(r'Median latitude (sea-ice min) [$^{\circ}$N]', fontsize=fontsize)
    plt.xlabel('Time', fontsize=fontsize)
    plt.title(titlestr, fontsize=fontsize)
    plt.savefig(figstr)
    plt.close(fig)

# -----------------------------------------------------------------------------
print('** END')
