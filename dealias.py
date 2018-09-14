import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import numpy as np
from pyart.io.nexrad_archive import read_nexrad_archive
import os

# fit cosine
def fit_cosine(xdata, ydata, freq):
    cproj = np.inner(np.cos(xdata*freq), ydata)
    sproj = np.inner(np.sin(xdata*freq), ydata)
    phi = np.arctan2(sproj, cproj)
    amp = np.sqrt(cproj**2.+sproj**2.)/len(xdata)
    yfit = amp*np.cos(xdata+phi)
    print cproj, sproj
    return yfit

# unfold segment
def unfold_seg(vel1d, azi1d, azi1, azi2, nyqv, sgn_fac):
    numazi = len(vel1d)
    print azi1, azi2, azi1d.shape, numazi
    ufld_ind = np.arange(numazi)[(azi1d>=azi1)&(azi1d<azi2)]
    vel_ufld = np.empty([numazi])
    vel_ufld[:] = vel1d
    print sgn_fac
    vel_ufld[ufld_ind] = sgn_fac*2.*nyqv+vel1d[ufld_ind]
    return vel_ufld

# open file
#-----------------------------------
print 'Opening file...'

yyyy = '2018'
mm = '09'
dd = '14'

hh = '12'
mn = '03'
ss = '46'
site = 'KLTX'

direc = 'aws_nexrad/'
filename = '{}{}{}{}_{}{}{}_V06'.format(site, yyyy, mm, dd,
                                        hh, mn, ss)
rad = read_nexrad_archive(direc+filename)
rad_sw = rad.extract_sweeps([1])

# get variables
#-----------------------------------
elev_p = rad_sw.elevation['data']
azi_p = 90.-rad_sw.azimuth['data']
ran = rad_sw.range['data']
vel_p = rad_sw.fields['velocity']['data']
nyvel = rad.get_nyquist_vel(1)

azi_p[azi_p<0.] = azi_p[azi_p<0.]+360.

dims = vel_p.shape
numradials = dims[0]+1
numgates = dims[1]

# expand radially to remove no data spike
elev = np.ma.empty([numradials])
azi = np.ma.empty([numradials])
vel = np.ma.empty([numradials, numgates])

elev[0:numradials-1] = elev_p
elev[numradials-1] = elev_p[0]
azi[0:numradials-1] = azi_p
azi[numradials-1] = azi_p[0]
vel[0:numradials-1,:] = vel_p
vel[numradials-1,:] = vel_p[0]

angle = np.mean(elev)

# mask velocity
#-----------------------------------------------
vel = np.ma.masked_where(vel==0., vel)
vel = np.ma.masked_where(vel==1., vel)

# calculate x and y coordinates (wrt beampath) for plotting
#-----------------------------------------------------------
ran_2d = np.tile(ran,(numradials,1))
azi.shape = (azi.shape[0], 1)
azi_2d = np.tile(azi,(1,numgates))

radz = 10.
erad = np.pi*angle/180.

ke = 4./3.
a = 6378137.

# beam height and beam distance
zcor = np.sqrt(ran_2d**2.+(ke*a)**2.+2.*ran_2d*ke*a*np.sin(erad))-ke*a+radz
scor = ke*a*np.arcsin(ran_2d*np.cos(erad)/(ke*a+zcor))/1000.

xcor = scor*np.cos(np.pi*azi_2d/180.)
ycor = scor*np.sin(np.pi*azi_2d/180.)

gtind = 243

# calculate change in vr between consecutive radials
chvr = vel[1:,gtind]-vel[0:-1,gtind]
chvr_sgn = chvr/np.abs(chvr)
indicator = np.zeros([numradials-1])
indicator[np.abs(chvr)>1.2*nyvel] = 1.*chvr_sgn[np.abs(chvr)>0.8*nyvel]
azi_seg = azi[:-1,0][indicator!=0]
ind_seg = indicator[indicator!=0]

# extend indicator and azimuth arrays at each end
numind = len(ind_seg)
ind_ex = np.zeros([numind+2])
azi_ex = np.zeros([numind+2])
ind_ex[1:-1] = ind_seg
ind_ex[0] = ind_seg[-1]
ind_ex[-1] = ind_seg[0]

azi_ex[1:-1] = azi_seg
azi_ex[0] = azi_seg[-1]
azi_ex[-1] = azi_seg[0]

# try unfolding one segment
print ind_ex
print azi_ex
consec_inds = []
for i in range(numind):
    ind_m1 = ind_ex[i]
    ind_z = ind_ex[i+1]
    ind_p1 = ind_ex[i+2]
    if((ind_m1==ind_z)|(ind_p1==ind_z)):
        consec_inds.append(i+1)

print consec_inds
# mark out unfolded areas and folded areas
numpair = len(consec_inds)/2
unfolded = np.zeros([numradials])-1.
for i in range(numpair):
    azi1 = azi_ex[consec_inds[i*2+1]]
    azi2 = azi_ex[consec_inds[i*2]]
    unfolded[(azi[:,0]>=azi1)&(azi[:,0]<azi2)] = 1.

# try unfolding
vel_org = vel[:,gtind].copy()
vel_test = np.empty([numradials])
vel_test[:] = vel[:,gtind]
sgn_fac = -vel_org[unfolded<0.]/np.abs(vel_org[unfolded<0.])
vel_test[unfolded<0.] = 2.*nyvel*sgn_fac+vel_org[unfolded<0.]

# fit cosine
freq = 1.
cfit = fit_cosine(azi[:,0]*np.pi/180., vel_test, freq)

'''
cind_str = consec_inds[0]
sgn_fac = 1.
for i in range(cind_str-1):
    azi1 = azi_ex[cind_str-i]
    azi2 = azi_ex[cind_str-i-1]
    vel[:,gtind] = unfold_seg(vel[:,gtind], azi[:,0], azi1, azi2, nyvel, sgn_fac)
    sgn_fac = -1*sgn_fac
'''

# plot
'''
plt.figure(0)
plt.pcolormesh(xcor, ycor, vel, cmap='spectral', vmin=-32., vmax=32.)
plt.plot(xcor[:,gtind], ycor[:,gtind], 'k--', lw=3.)
ax = plt.gca()
ax.set_aspect(1.)
ax.set_xlim([-250., 250.])
ax.set_ylim([-250., 250.])
plt.savefig('velocity.png')
'''
plt.figure(1)
plt.plot(azi, vel_test, 'b--', lw=3.)
plt.plot(azi, cfit, 'g--', lw=3.)
plt.plot(azi, vel_org, 'k--', lw=3.)
plt.plot(azi[1:], indicator*10, 'r--', lw=3.)
plt.scatter(azi, vel[:,gtind], c=unfolded, cmap='spectral', s=40.)
plt.colorbar()
plt.savefig('circle.png')
