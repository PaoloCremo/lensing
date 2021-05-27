#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 11:36:31 2021

@author: paolo
"""
#%%load packages and functions
# from Dropbox.PhD.Python.program_py.packages.pyCBC_function import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#from Dropbox.PhD.Python.program_py.packages import LISA as li
#lisa = li.LISA()
import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
#%%
# define constants
H0 = 74.
Om = 0.3061
Og = 0 #10**(-5)
Ok = 0
Ol = 1. - Om - Og - Ok
w_0 = -1.
w_a = 0.
c = 299_792_459 #m/s
G = 6.67408*10**(-11) #m^3 kg^-1 s^-2
smtokg=1.98847*10**30
pctomt = 3.08567758 * 10**16 #pc --> meters
TSUN    = 4.92549232189886339689643862e-6 # mass of sun in seconds (G=C=1)
YEAR   = 3.15581497632e7  # year in seconds

# define param for plot
fontSz = 15
fontsz = 13
fontssz = 11

df_clr = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']



#%% define models
models = ['SIS', 'NFW', 'NFW_2', 'GNFW_2']
def rho_SIS(r, sigma):
    rr = r**2*sigma**2/(2*np.pi*G*r**2)

    return rr

def rho_NFW(r, rhos, rs, gamma):
    rr = r**2*rhos/(((r/rs)**gamma)*(1+(r/rs))**(3-gamma))
    
    return rr

#%% parameters

sigma = 66_260

rhos_NFW = 7.95*10**-22
rs_NFW = 29.38*pctomt*10**3

rhos_NFW2 = rhos_NFW
rs_NFW2 = 11.12*pctomt*10**3

rhos_GNFW2 = 4.94*10**-22
rs_GNFW2 = 4.92*pctomt*10**3

#%%
rs = np.concatenate((np.arange(10**17, 10**19, 10**17), np.arange(10**19, 10**21, 10**19), np.arange(10**21, 10**23, 10**21)))
print(len(rs))

y_SIS = [4*np.pi*integrate.quad(rho_SIS, 0, r, args=(sigma))[0]/smtokg for r in rs]
y_NFW = [4*np.pi*integrate.quad(rho_NFW, 0, r, args=(rhos_NFW, rs_NFW, 1))[0]/smtokg for r in rs]
y_NFW2 = [4*np.pi*integrate.quad(rho_NFW, 0, r, args=(rhos_NFW2, rs_NFW2, 1))[0]/smtokg for r in rs]
y_GNFW2 = [4*np.pi*integrate.quad(rho_NFW, 0, r, args=(rhos_GNFW2, rs_GNFW2, 2))[0]/smtokg for r in rs]


plt.figure()
plt.plot(rs/pctomt, y_SIS, c='r', label=models[0])
plt.plot(rs/pctomt, y_NFW, c='b', label=models[1])
plt.plot(rs/pctomt, y_NFW2, c=df_clr[9], label='NFW-2')
plt.plot(rs/pctomt, y_GNFW2, c='green', label='gNFW$_\gamma=2$')

plt.xlabel('r [pc]', fontsize=fontsz)
plt.ylabel('M [$M_\odot$]', fontsize=fontsz)

plt.xscale('log')
plt.yscale('log')

plt.legend(ncol=4, loc='lower center', bbox_to_anchor=(0.5,1.), fontsize=fontsz,framealpha=0.5 )
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)
#plt.legend()
# plt.show()
#%%
# '''
#plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/models_masses.png',
plt.savefig('/Users/paolocremonese/Dropbox/PhD/plots/lens_masses/models_masses.png',
            dpi=300, format='png', transparent=True, bbox_inches="tight")
# plt.close()
# '''
