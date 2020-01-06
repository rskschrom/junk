from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from siphon.catalog import TDSCatalog
import numpy as np
import os
import datetime
import wind as wnd

# interpolate data to height agl
def interpolate_agl(data_iso, hgt_levs, hgt_agl_iso):
    nhgt = len(hgt_levs)
    dims = data_iso.shape
    data_agl = np.ma.masked_all([dims[0],nhgt,dims[2],dims[3]])
    dh = hgt_levs[1]-hgt_levs[0]
    for i in range(nhgt):
        wgt = np.exp(-(hgt_agl_iso-hgt_levs[i])**2./(dh)**2.)
        wgt = wgt/np.sum(wgt, axis=1)
        data_agl[:,i,:,:] = np.sum(data_iso*wgt, axis=1)

    return data_agl

# read in latest gfs data
best_hrrr = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/'+
                      'grib/NCEP/RAP/CONUS_20km/latest.html')
#best_hrrr = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/'+
#                      'grib/NCEP/GFS/Global_0p5deg/latest.html')
ref_time = best_hrrr.metadata['documentation']['Reference Time'][0]
yyyy = ref_time[0:4]
mm = ref_time[5:7]
dd = ref_time[8:10]
hh = ref_time[11:13]
print(yyyy, mm, dd, hh)

time = datetime.datetime(year=int(yyyy), month=int(mm), day=int(dd), hour=int(hh))
time.replace(tzinfo=datetime.timezone.utc)
best_ds = best_hrrr.datasets.filter_time_nearest(time, regex=None)
ncss = best_ds.subset()
query = ncss.query()
#print(ncss.variables)

# want time
fhour = 2
forecast_time = time+datetime.timedelta(hours=fhour)

# get subset
print('downloading data...')
query.lonlat_box(-120.,-60.,25.,55.)
#query.lonlat_box(-95.,-65.,32.,50.)
query.variables('Geopotential_height_surface',
                'Geopotential_height_isobaric',
                'u-component_of_wind_isobaric',
                'v-component_of_wind_isobaric',
                'Temperature_isobaric').accept('netcdf')
query.add_lonlat(True)
#query.all_times()
query.time(forecast_time)
data = ncss.get_data(query)

# get variable fields
geo_iso = data.variables['Geopotential_height_isobaric'][:]
geo_sfc = data.variables['Geopotential_height_surface'][:]
u_iso = data.variables['u-component_of_wind_isobaric'][:]
v_iso = data.variables['v-component_of_wind_isobaric'][:]
tmp_iso = data.variables['Temperature_isobaric'][:]
plev = data.variables['isobaric'][:]
lat = data.variables['lat'][:]
lon = data.variables['lon'][:]

# get data on height levels
geo_agl = geo_iso-geo_sfc
dh = 125
hgt_levs = np.arange(0., 7000.+dh, dh)
print(hgt_levs)

#print(hgt_levs)
tmp_agl = interpolate_agl(tmp_iso, hgt_levs, geo_agl)
u_agl = interpolate_agl(u_iso, hgt_levs, geo_agl)
v_agl = interpolate_agl(v_iso, hgt_levs, geo_agl)

plev.shape = (1,len(plev),1,1)
pres = np.tile(plev, (1,1,geo_iso.shape[2],geo_iso.shape[3]))
pres_agl = interpolate_agl(pres, hgt_levs, geo_agl)
dens_agl = pres_agl/tmp_agl
wgt = pres_agl/np.sum(pres_agl, axis=1)
print(pres_agl[0,:,67,123])

# calculate storm-relative helicity
u_sm, v_sm, u_sh, v_sh = wnd.storm_motion(u_agl, v_agl, hgt_levs, wgt=wgt)
srh = wnd.storm_relative_helicity(u_agl, v_agl, hgt_levs, wgt=wgt)

# map plot
fig = plt.figure(figsize=(12,10))
ax = fig.add_subplot(1,1,1, projection=ccrs.LambertConformal())

# plot filled contours
skf = 5
srh_levs = 50*np.arange(10)+50
cols = []
for slev in srh_levs:
    if slev<200.:
        cols.append('lightskyblue')
    elif (slev>=200.)&(slev<300.):
        cols.append('dodgerblue')
    else:
        cols.append('navy')

# label contours
cs = ax.contour(lon, lat, srh, levels=srh_levs, colors=cols, linewidths=1., transform=ccrs.PlateCarree())
cls = ax.clabel(cs, inline=True, inline_spacing=2., fmt='%3.0f', fontsize=12)
for cl in cls:
    cl.set_rotation(0)

ax.barbs(lon[::skf,::skf], lat[::skf,::skf], 1.94*u_sm[::skf,::skf], 1.94*v_sm[::skf,::skf],
         color='goldenrod', length=4., linewidths=0.5, transform=ccrs.PlateCarree())

# map features
maxlon = -90.
minlon = -70.
minlat = 35.
maxlat = 45.
ax.set_extent([minlon,maxlon,minlat,maxlat], crs=ccrs.PlateCarree())
shapename = 'admin_1_states_provinces_lakes_shp'
states_shp = shpreader.natural_earth(resolution='50m',
                                     category='cultural', name=shapename)
ax.add_geometries(shpreader.Reader(states_shp).geometries(), ccrs.PlateCarree(), edgecolor='gray',
                  facecolor='none')

# set axis titles
ax.set_title(f'Forecast time: {forecast_time} UTC',
             x=0.0, y=1.02, horizontalalignment='left', fontsize=16)
plt.suptitle(f'{hh}-UTC HRRR forecast - {yyyy}-{mm}-{dd}', y=0.98, fontsize=26)
imgname = 'hrrr_srh.png'
plt.savefig(imgname)

