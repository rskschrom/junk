import numpy as np
from scipy import integrate

# integrate ellipsoid shape factors numerically
def ellipsoid_shape_facs(a, b, c):
    f = lambda x,av,bv,cv,dv : 0.5*av*bv*cv/((x+dv**2.)*
                               np.sqrt((x+av**2.)*(x+bv**2.)*(x+cv**2.)))

    la, err = integrate.quad(f, 0., np.inf, args=(a,b,c,a))
    lb, err = integrate.quad(f, 0., np.inf, args=(a,b,c,b))
    lc, err = integrate.quad(f, 0., np.inf, args=(a,b,c,c))

    return la, lb, lc

# ellipsoid polarizabilities
def ellipsoid_polz(diel, a, b, c):
    la, lb, lc = ellipsoid_shape_facs(a, b, c)
    alph_a = 4./3.*np.pi*a*b*c*(diel-1.)/(1.+la*(diel-1.))
    alph_b = 4./3.*np.pi*a*b*c*(diel-1.)/(1.+lb*(diel-1.))
    alph_c = 4./3.*np.pi*a*b*c*(diel-1.)/(1.+lc*(diel-1.))

    return alph_a, alph_b, alph_c

# point spin orientation
def point_spin(phi_p, theta_p, phi_s, alpha_a, alpha_b, alpha_c):
    # create point rotation matrix
    cpp = np.cos(phi_p)
    ctp = np.cos(theta_p)
    spp = np.sin(phi_p)
    stp = np.sin(theta_p)
    apoint = np.array([[cpp*ctp,spp*ctp,-stp],
                       [-spp,cpp,0.],
                       [cpp*stp,spp*stp,ctp]])

    # create spin rotation matrix
    cps = np.cos(phi_s)
    sps = np.sin(phi_s)
    aspin = np.array([[cps,-sps,0.],
                       [sps,cps,0.],
                       [0.,0.,1.]])

    # transform polarizability tensor
    alpha_tensor = np.array([[alpha_a,0.,0.],
                             [0.,alpha_b,0.],
                             [0.,0.,alpha_c]])
    arot = np.matmul(aspin, apoint)
    alpha_rot = np.matmul(arot.T, np.matmul(alpha_tensor, arot))
    return alpha_rot

# scattering matrix from incident and scattering direction
# th_i - incident elevation angle from x-y plane
# ph_i - incident azimuthal angle from x axis
# th_s - scattering elevation angle from x-y plane
# ph_s - scattering azimuthal angle from x axis
def scat_matrix(th_i, phi_i, th_s, phi_s, alpha_tensor, k):
    sti = np.sin(th_i)
    cti = np.cos(th_i)
    spi = np.sin(phi_i)
    cpi = np.cos(phi_i)
    apol = np.array([[-spi, cpi, 0.],
                     [-cpi*sti, -spi*sti, cti]])
    alpha_pol = np.matmul(np.matmul(apol, alpha_tensor), apol.T)
    smat = k**2./(4.*np.pi)*alpha_pol
    return smat

# get shape factors for an ellipsoid
a = 2.
b = 1.
c = 0.1
diel = complex(3.17, 0.025)
wavl = 32.1
k = 2.*np.pi/wavl
alpa, alpb, alpc = ellipsoid_polz(diel, a, b, c)

# try rotations and get scattering matrix
phi_p = 0.
theta_p = 0.
phi_s = 0.
alp_rot = point_spin(phi_p, theta_p, phi_s, alpa, alpb, alpc)

el_i = 0.
az_i = 0.
el_s = 0.
az_s = 0.
smat = scat_matrix(el_i, az_i, el_s, az_s, alp_rot, k)
shh = smat[0,0]
shv = smat[0,1]
svv = smat[1,1]
zdr = 10.*np.log10(np.abs(shh)**2./np.abs(svv)**2.)
print zdr
