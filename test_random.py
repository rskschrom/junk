import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import scipy.stats as st
from scipy.special import erfc

# get frequency histogram
def freq_hist(minb, maxb, nbin, data):
    bins = np.linspace(minb, maxb, nbin+1)
    freq = np.empty([nbin])

    for i in range(nbin):
        bin_data = data[(data>bins[i])&(data<=bins[i+1])]
        freq[i] = len(bin_data)/float(ndata)

    return bins, freq

# open file
ncfile = Dataset('nsakazrgeC1.b1.20130502.000003.cdf', 'r')
print ncfile
rxnoise = ncfile.variables['rx_noise'][:]
zh = ncfile.variables['reflectivity_copol'][:]

# get random noise
noise = zh[:,550]
noise = (noise-np.mean(noise))/np.std(noise)
ndata = len(noise)

# bin data
nbin = 20
minb = -4.5
maxb = 4.5
bins, freq = freq_hist(minb, maxb, nbin,  noise)

# get idealized values
f_ideal = np.empty([nbin])
for i in range(nbin):
    f_ideal[i] = st.norm.cdf(bins[i+1])-st.norm.cdf(bins[i])

# do chi^2 test
chi2, p = st.chisquare(freq, f_ideal)
print p, chi2

# do the same with uniform transformation
u = 0.5*erfc(-noise/np.sqrt(2.))
nbin = 20
minb = 0.
maxb = 1.
bins, freq = freq_hist(minb, maxb, nbin,  u)
print freq
chi2, p = st.chisquare(freq, np.full([nbin], 1./nbin))
print p, chi2

# plot
plt.plot(noise, 'k-', lw=1.)
plt.savefig('rx.png')
