import pygrib
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# open file
grbs = pygrib.open('hrrr.t23z.wrfsfcf17.grib2')
grbs.seek(0)
#for grb in grbs:
#    print grb

geolist = grbs.select(name='Precipitable water')

# get lat and lon
lat2d, lon2d = geolist[0].latlons()
nlat = lat2d.shape[0]
nlon = lat2d.shape[1]

# create array
geo = geolist[0]
print geolist
field = geo.values
print field.shape

# plot
#map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
map = Basemap(projection='merc',llcrnrlat=20.,urcrnrlat=55.,\
              llcrnrlon=-125.,urcrnrlon=-60.,lat_ts=20,resolution='l')
map.drawcoastlines(linewidth=0.5)
map.drawcountries(linewidth=0.5)
map.drawcounties(linewidth=0.2)
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))

x, y = map(lon2d, lat2d)
map.contourf(x, y, field, cmap='viridis')
plt.colorbar()
plt.savefig('test.png')
