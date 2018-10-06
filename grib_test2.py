import pygrib
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

# open file
grbs = pygrib.open('test.grib')
grbs.seek(0)
geolist = grbs.select(name='Geopotential')
nmem = 51
ntime = len(geolist)/nmem

# get lat and lon
lat2d, lon2d = geolist[0].latlons()
nlat = lat2d.shape[0]
nlon = lat2d.shape[1]

# create arrays of geopotential over time and ensemble member
geo = np.empty([nmem,ntime,nlat,nlon])
for i in range(nmem):
    print i
    for j in range(ntime):
        geo[i,j,:,:] = geolist[i*ntime+j].values

print geo.shape

# plot geopotential of first member at first time
map = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
map.drawcoastlines(linewidth=0.5)
map.drawcountries(linewidth=0.5)
map.drawmeridians(np.arange(0,360,30))
map.drawparallels(np.arange(-90,90,30))

x, y = map(lon2d, lat2d)
map.contourf(x, y, geo[0,0,:,:], cmap='Spectral')
plt.colorbar()
plt.savefig('test.png')

