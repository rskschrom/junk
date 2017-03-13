from pydap.client import open_url
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

dataset = open_url('http://nomads.ncep.noaa.gov:9090/dods/gfs_1p00/gfs20170313/gfs_1p00_06z')
print dataset.keys()
pcp = dataset['apcpsfc']
acum_precip = pcp['apcpsfc']
lat = pcp['lat']
lon = pcp['lon']
end_precip = acum_precip[-1,:,:].transpose().squeeze()

lon2d, lat2d = np.meshgrid(lon, lat, indexing='ij')

print lon2d.shape, lat2d.shape, end_precip.shape
loni1 = 240
loni2 = 300
lati1 = 90
lati2 = 150

print lat[lati1], lat[lati2]
print lon[loni1], lon[loni2]

# plot
plt.contourf(lon2d[loni1:loni2,lati1:lati2],
             lat2d[loni1:loni2,lati1:lati2],
             end_precip[loni1:loni2,lati1:lati2])
plt.savefig('opendap_test.png')

