from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from siphon.catalog import TDSCatalog
import numpy as np
import os

# read in latest gfs data
best_hrrr = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/'+
                      'grib/NCEP/HRRR/CONUS_2p5km/latest.html')
ref_time = best_hrrr.metadata['documentation']['Reference Time'][0]
yyyy = ref_time[0:4]
mm = ref_time[5:7]
dd = ref_time[8:10]
hh = ref_time[11:13]
print yyyy, mm, dd, hh
best_ds = best_hrrr.datasets.itervalues().next()
ncss = best_ds.subset()
query = ncss.query()

# get subset
print 'downloading data...'
now = datetime.utcnow()
query.lonlat_box(-90.,-70.,40.,45.).time(now+timedelta(hours=0))
query.variables('Reflectivity_height_above_ground').accept('netcdf')
query.add_lonlat(True)
#query.vertical_level(50000.)
query.all_times()
data = ncss.get_data(query)

# get variable fields
ref = data.variables['Reflectivity_height_above_ground']
lat = data.variables['lat']
lon = data.variables['lon']

# get dimensions and forecast hours
dims = ref.shape
print dims
nfh = dims[0]

# plot stuff out of loop
rmi = -10.
rma = 80.
ntick = 19.
step = (rma-rmi)/(ntick-1)
ticks = np.linspace(rmi, rma+2*step, ntick+2)

levs = np.linspace(rmi, rma, ntick)

# loop over forecast hours
for fh in range(nfh):

    # plot
    print 'plotting fh {:02d}'.format(fh)
    fig = plt.figure(fh)
    ax = fig.add_subplot(1,1,1)
    m = Basemap(width=1400000,height=1100000,rsphere=6371229.0,
                projection='lcc',lat_1=25.,lat_2=25.,lon_0=284.,lat_0=41.,
                resolution ='h',area_thresh=1000.)

    x, y = m(lon[:], lat[:])

    # plot filled contours
    cf = m.contourf(x, y, ref[fh,0,:,:], cmap='Spectral_r', levels=levs, lw=0., vmin=rmi, vmax=rma)
    cbar = m.colorbar(cf, pad='7%')
    cb_la = ['{:.1f}'.format(ti) for ti in ticks]
    cbar.set_ticks(ticks)
    cbar.ax.set_yticklabels(cb_la, fontsize=18)
    cbar.set_label('dBZ', fontsize=22)

    # draw coastlines, meridians and parallels.
    m.drawcoastlines(linewidth=1.5, color='k')
    m.drawcountries(linewidth=1., color='k')
    m.drawcounties(linewidth=0.6, color='k')
    m.drawstates(linewidth=1., color='k')

    # set axis titles
    ax.set_title('Forecast hour: {:02d} - 1-km REF (shaded)'.format(fh),
                 x=0.0, y=1.02, horizontalalignment='left', fontsize=20)
    plt.suptitle('{}-UTC HRRR forecast - {}{}{}'.format(hh, yyyy, mm, dd), y=0.92, fontsize=26)
    #plt.tight_layout()
    imgname = 'ref_hrrr_{:02d}.png'.format(fh)
    plt.savefig(imgname)
    os.system('convert -trim {} {}'.format(imgname, imgname))
    plt.close()
