from pybufrkit.decoder import Decoder
from pybufrkit.renderer import FlatJsonRenderer
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cm
from matplotlib import patches
from mpl_toolkits.basemap import Basemap
import os

# colorblind radar map from pyart
#----------------------------------
def rgb_rainbow_24(nc):
    path1 = np.linspace(0.8*np.pi, 1.8*np.pi, nc)
    path2 = np.linspace(-0.33*np.pi, 0.33*np.pi, nc)

    y = np.concatenate([np.linspace(0.3, 0.85, nc*2//5),
                        np.linspace(0.9, 0.0, nc - nc*2//5)])
    u = 0.40*np.sin(path1)
    v = 0.55*np.sin(path2) + 0.1

    rgb_from_yuv = np.array([[1, 0, 1.13983],
                             [1, -0.39465, -0.58060],
                             [1, 2.03211, 0]])
    cdata = np.empty([nc,3])
    for i in range(len(y)):
        r = 1.*y[i]+0.*u[i]+1.13983*v[i]
        g = 1.*y[i]-0.39465*u[i]-0.5806*v[i]
        b = 1.*y[i]+2.03211*u[i]+0.*v[i]
        cdata[i,0] = min(1.,r)
        cdata[i,1] = min(1.,g)
        cdata[i,2] = min(1.,b)

    cmap = cm.ListedColormap(np.abs(cdata), 'Hoymeyer_cb')
    return cmap


# convert list of ensemble data to dictionary
def ens_data_dict(ens_tc_arr):
    nmem = len(ens_tc_arr)
    ntime = 41   
    ftime = np.arange(ntime)*6
    lat = np.empty([ntime,nmem])
    lon = np.empty([ntime,nmem]) 
    msl = np.empty([ntime,nmem])
    wspd = np.empty([ntime,nmem])

    # loop through ensemble members
    for i in range(nmem):
        header = ens_tc_arr[i][0:25]
        lat[0,i] = header[17]
        lon[0,i] = header[18]
        msl[0,i] = header[19]
        wspd[0,i] = header[21]
        #print lat[0,i], lon[0,i], msl[0,i], wspd[0,i]

        forecast = ens_tc_arr[i][25:]
        lat[1:,i] = np.array(forecast[3::10])
        lon[1:,i] = np.array(forecast[4::10])
        msl[1:,i] = np.array(forecast[5::10])
        wspd[1:,i] = np.array(forecast[7::10])

    # assemble dictionary
    lat = np.ma.masked_invalid(lat)
    lon = np.ma.masked_invalid(lon)
    msl = np.ma.masked_invalid(msl)
    wspd = np.ma.masked_invalid(wspd)

    print header
    name = header[4]
    year = header[8]
    month = header[9]
    day = header[10]
    hour = header[11]
    cen_lat = header[14]
    cen_lon = header[15]
    cen_msl = header[19]
    cen_wspd = header[21]

    ens_tc_dict = {'storm_metadata':{'name':name, 'year':year, 'month':month,
                                     'day':day, 'hour':hour, 'initial_lat':cen_lat,
                                     'initial_lon':cen_lon, 'init_msl':cen_msl,
                                     'init_wspd':cen_wspd},
                   'forecast_time':ftime, 'forecast_lat':lat,
                   'forecast_lon':lon, 'forecast_msl':msl, 'forecast_wspd':wspd}
    return ens_tc_dict

# read in bufr file
dec = Decoder()
dd = '07'
mm = '09'
yy = '18'
hh = '12'
name = 'FLORENCE'
latf = '-51p4degW'
lonf = '25degN'
s1 = 'A_JSXX01ECEP{}{}00'.format(dd, hh)
s2 = '_C_ECMP_20{}{}{}{}0000_tropical_cyclone_track_'.format(yy, mm, dd, hh)
s3 = '{}_{}_{}_bufr4.bin'.format(name, latf, lonf)
fn = s1+s2+s3

# read in file
print 'opening file...'
with open(fn, 'rb') as ins:
    bufr_message = dec.process(ins.read())

# convert to json
print 'conver to json...'
json_data = FlatJsonRenderer().render(bufr_message)
nj = len(json_data)
jd = [None]*nj
for i, j in enumerate(json_data):
    jd[i] = j

# put ensemble data into a dictionary
ens_data = jd[4][2]
nens = len(ens_data)
ens_tc_arr = [None]*nens

print 'making dictionary...'
for i in range(nens):
    ens_tc_arr[i] = ens_data[i]

ens_dict = ens_data_dict(ens_tc_arr)
lat = ens_dict['forecast_lat']
lon = ens_dict['forecast_lon']
wspd = ens_dict['forecast_wspd']*1.944
msl = ens_dict['forecast_msl']/1.e2
time = ens_dict['forecast_time']
ntime = len(time)

# plot
fig = plt.figure(0)
ax = fig.add_subplot(1,1,1)
m = Basemap(llcrnrlon=-100.,llcrnrlat=18.,urcrnrlon=-40.,urcrnrlat=52.,
            projection='lcc',lat_1=10.,lat_2=30.,lon_0=-75.,
            resolution ='l',area_thresh=1000.)

ncol = 31
minc = 920.
maxc = 1010
cmap = rgb_rainbow_24(ncol)
colors = np.empty([ntime,nens,4])

# normalize colors by msl
val_norm = (msl-minc)/(maxc-minc)
for i in range(ntime):
    for j in range(nens):
        colors[i,j,:] = cmap(val_norm[i,j])

for j in range(nens):
    x, y = m(lon[:,j], lat[:,j])
    for i in range(ntime-1):
        xseg = x[i:i+2]
        yseg = y[i:i+2]
        m.plot(xseg, yseg, color=colors[i+1,j,:], lw=1., alpha=0.5)


# plot track locations
x, y = m(lon, lat)
msp = m.scatter(x, y, c=msl, s=30., lw=0., cmap=cmap, vmin=minc, vmax=maxc)

# dummy colormap
cb = m.colorbar(msp, location='right', pad='10%')
cb.set_label('MSLP (hPa)', fontsize=22)

# calculate covariance matrices and ellipses at each time
for i in range(ntime):
    xe, ye = m(lon[i,:], lat[i,:])
    cov_mat = np.ma.cov(xe, ye)

    eigv, ueigv = np.linalg.eig(cov_mat)
    vx = eigv[0]*ueigv[:,0]
    vy = eigv[1]*ueigv[:,1]
    
    ang = np.arctan2(vx[1], vx[0])
    print ang*180./np.pi

    xc = np.ma.mean(xe)
    yc = np.ma.mean(ye)

    e1 = patches.Ellipse((xc, yc), np.sqrt(eigv[0])*4., np.sqrt(eigv[1])*4., angle=ang*180./np.pi,
                         edgecolor='k', linewidth=4., fill=False, facecolor='grey', alpha=0.3)
    print 'adding patch', i
    ax.add_patch(e1)

# plot median track
latmed = np.ma.median(lat, axis=1)
lonmed = np.ma.median(lon, axis=1)
xm, ym = m(lonmed, latmed)
m.plot(xm, ym, 'k-', lw=4.)

# draw coastlines, meridians and parallels.
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1.)
m.drawstates(linewidth=1.)
m.drawparallels(np.arange(10,70,20),labels=[1,1,0,0])
m.drawmeridians(np.arange(-100,0,20),labels=[0,0,0,1])

# set axis titles
ax.set_title('{}-UTC ECMWF Ensemble - {} September 20{}'.format(hh, dd, yy),
             x=0.38, y=1.02, fontsize=26)

plt.tight_layout()
imgname = 'ellipses.png'
plt.savefig(imgname)
os.system('convert -trim {} {}'.format(imgname, imgname))
