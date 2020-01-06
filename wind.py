import numpy as np

# calculate storm motion (Bunkers et al. 2000)
def storm_motion(u_agl, v_agl, hgt_levs, wgt=None):
    # get indicies for height level (should have every 500m agl)
    ind00 = np.argmin(np.abs(hgt_levs-0.))
    ind05 = np.argmin(np.abs(hgt_levs-500.))+1
    ind55 = np.argmin(np.abs(hgt_levs-5500.))
    ind60 = np.argmin(np.abs(hgt_levs-6000.))+1

    # calculate shear
    umn_00_05 = np.mean(u_agl[0,ind00:ind05,:,:], axis=0)
    umn_55_60 = np.mean(u_agl[0,ind55:ind60,:,:], axis=0)
    vmn_00_05 = np.mean(v_agl[0,ind00:ind05,:,:], axis=0)
    vmn_55_60 = np.mean(v_agl[0,ind55:ind60,:,:], axis=0)
    u_shear = umn_55_60-umn_00_05
    v_shear = vmn_55_60-vmn_00_05
    mag_shear = np.sqrt(u_shear**2.+v_shear**2.)
    ush_dir = u_shear/mag_shear
    vsh_dir = v_shear/mag_shear

    # calculate mean wind from 0-6 km
    u_mn = np.mean(u_agl[0,ind00:ind60,:,:], axis=0)
    v_mn = np.mean(v_agl[0,ind00:ind60,:,:], axis=0)

    if wgt is not None:
        u_mn = np.sum(wgt[0,ind00:ind60,:,:]*u_agl[0,ind00:ind60,:,:], axis=0)
        v_mn = np.sum(wgt[0,ind00:ind60,:,:]*v_agl[0,ind00:ind60,:,:], axis=0)

    # deviation from mean wind
    u_shft = vsh_dir*7.5
    v_shft = -ush_dir*7.5
    u_sm = u_mn+u_shft
    v_sm = v_mn+v_shft
    #u_sm = u_mn[:]
    #v_sm = v_mn[:]
    return u_sm, v_sm, u_shear, v_shear

# calculate storm motion for 1d data (Bunkers et al. 2000)
def storm_motion1d(u_agl, v_agl, hgt_levs, wgt=None):
    # get indicies for height level (should have every 500m agl)
    ind00 = np.argmin(np.abs(hgt_levs-0.))
    ind05 = np.argmin(np.abs(hgt_levs-500.))+1
    ind55 = np.argmin(np.abs(hgt_levs-5500.))
    ind60 = np.argmin(np.abs(hgt_levs-6000.))+1

    # calculate shear
    umn_00_05 = np.mean(u_agl[ind00:ind05])
    umn_55_60 = np.mean(u_agl[ind55:ind60])
    vmn_00_05 = np.mean(v_agl[ind00:ind05])
    vmn_55_60 = np.mean(v_agl[ind55:ind60])
    u_shear = umn_55_60-umn_00_05
    v_shear = vmn_55_60-vmn_00_05
    mag_shear = np.sqrt(u_shear**2.+v_shear**2.)
    ush_dir = u_shear/mag_shear
    vsh_dir = v_shear/mag_shear

    # calculate mean wind from 0-6 km
    u_mn = np.mean(u_agl[ind00:ind60])
    v_mn = np.mean(v_agl[ind00:ind60])

    if wgt is not None:
        wgt = wgt[ind00:ind60]/np.sum(wgt[ind00:ind60])
        u_mn = np.sum(wgt[ind00:ind60]*u_agl[ind00:ind60])
        v_mn = np.sum(wgt[ind00:ind60]*v_agl[ind00:ind60])

    # deviation from mean wind
    u_shft = vsh_dir*7.5
    v_shft = -ush_dir*7.5
    u_sm = u_mn+u_shft
    v_sm = v_mn+v_shft
    return u_sm, v_sm

# calculate storm-relative helicity (M+R 2011)
def storm_relative_helicity(u_agl, v_agl, hgt_levs, srh_depth=3000., wgt=None):
    depth_ind = np.argmin(np.abs(hgt_levs-srh_depth))+1
    u_sm, v_sm,_,_ = storm_motion(u_agl, v_agl, hgt_levs, wgt)
    srh = np.sum((u_agl[0,1:depth_ind,:,:]-u_sm)*(v_agl[0,:depth_ind-1,:,:]-v_sm)-\
                 (u_agl[0,:depth_ind-1,:,:]-u_sm)*(v_agl[0,1:depth_ind,:,:]-v_sm), axis=0)
    return srh

# calculate storm-relative helicity for 1d data (M+R 2011)
def storm_relative_helicity1d(u_agl, v_agl, hgt_levs, srh_depth=3000., wgt=None):
    depth_ind = np.argmin(np.abs(hgt_levs-srh_depth))+1
    u_sm, v_sm = storm_motion1d(u_agl, v_agl, hgt_levs, wgt)
    srh = np.sum((u_agl[1:depth_ind]-u_sm)*(v_agl[:depth_ind-1]-v_sm)-\
                 (u_agl[:depth_ind-1]-u_sm)*(v_agl[1:depth_ind]-v_sm))
    return srh

# interpolate 1d data to height agl
def interpolate_agl1d(data_iso, hgt_levs, hgt_agl_iso):
    nhgt = len(hgt_levs)
    data_agl = np.ma.masked_all([nhgt])
    for i in range(nhgt):
        if (i<nhgt-1):
            dh = hgt_levs[i+1]-hgt_levs[i]
        else:
            dh = hgt_levs[i]-hgt_levs[i-1]
        wgt = np.exp(-(hgt_agl_iso-hgt_levs[i])**2./(dh*1.2)**2.)
        wgt = wgt/np.sum(wgt)
        data_agl[i] = np.sum(data_iso*wgt)

    return data_agl
