#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 13 11:36:31 2021

@author: paolo
"""
#%%load packages and functions
from Dropbox.PhD.Python.program_py.packages.pyCBC_function import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from Dropbox.PhD.Python.program_py.packages import LISA as li
lisa = li.LISA()
lbl = ['SIS', 'NFW', '$gNFW_{\gamma=2}$', 'NFW-2', 'NFW-2\n$z_s=0.5$\n$z_L=0.25$']
#%% define case
case = 'h'#'d'
df = 10**-6
dt = 1
zS = 1.
DS = cosmo.luminosity_distance(zS).value #distance of the source in Mpc
zL = 0.5
Ml = 10**9
f_min = 10**-5
MTOT = 10**8*(1+zS)
q = 1.
m1 = MTOT/(1+q)
m2 = MTOT-m1
#%% time domain
l = 3157032
cf = 0#np.pi/4

hp, hc = get_td_waveform(approximant='IMRPhenomPv3', mass1=m1, mass2=m2,
                         delta_t=1, f_lower=f_min, inclination=np.pi/3,
                         # s1z = 0.7, s2z = 0.2,
                         distance=DS, coa_phase=cf)
hp.resize(l)
# l = len(hp)
ht = make_transform_with_C_r2c(hp.data)
fr = np.arange(l/2+1)/(l*dt)
# '''
#% degenerate wf
# '''
zS_deg = 0.55
DS_deg = cosmo.luminosity_distance(zS_deg).value 
hp_deg, hc_deg = get_td_waveform(approximant='IMRPhenomPv3', mass1=m1, mass2=m2,
                         delta_t=1, f_lower=f_min, inclination=np.pi/3,
                          # s1z = 0.7, s2z = 0.2,
                         distance=DS_deg, coa_phase=0)
hp_deg.resize(l)
# l = len(hp)
ht_deg = make_transform_with_C_r2c(hp_deg.data)
#% degenerate wf 2
# '''
zS_deg2 = 0.8
DS_deg2 = cosmo.luminosity_distance(zS_deg2).value 
hp_deg2, hc_deg = get_td_waveform(approximant='IMRPhenomPv3', mass1=m1, mass2=m2,
                         delta_t=1, f_lower=f_min, inclination=np.pi/3,
                          # s1z = 0.7, s2z = 0.2,
                         distance=DS_deg2, coa_phase=0)
hp_deg2.resize(l)
# l = len(hp)
ht_deg2 = make_transform_with_C_r2c(hp_deg2.data)
#%% load amps
l1 = 'SIS'
y1 = '1'
M1 = '10-9'
zL1 = '05'

af1 = '/home/paolo/Desktop/waveform/case_h/'+l1+'/amps/amps_case_h_H74_y'+y1+'_M'+M1+'_zL'+zL1+'_norm'
amps_1 = np.fromfile(af1, np.complex64)
amps_1 = np.concatenate((amps_1, np.zeros(int(l/2+1)-len(amps_1))))

hLt1 = ht * np.conj(amps_1)
# hLt1 = ht_deg * np.conj(amps_1)
hL1 = make_transform_with_C_c2r(hLt1)

#%  load amps 2
l2 = 'NFW'
y2 = '001061'
# y2 = '002804'
M2 = '10-9'

af2 = '/home/paolo/Desktop/waveform/case_h/'+l2+'/amps/amps_case_h_H74_y'+y2+'_M'+M2+'_zL05_norm'
amps_2 = np.fromfile(af2, np.complex64)
amps_2 = np.concatenate((amps_2, np.zeros(int(l/2+1)-len(amps_2))))

hLt2 = ht * np.conj(amps_2)
# hLt2 = ht_deg * np.conj(amps_2)
hL2 = make_transform_with_C_c2r(hLt2)
#%  load amps 3
l3 = 'GNFW_2'
y3 = '00634'
M3 = '10-9'

af3 = '/home/paolo/Desktop/waveform/case_h/'+l3+'/amps/amps_case_h_H74_y'+y3+'_M'+M3+'_zL05'#_FP'
amps_3 = np.fromfile(af3, np.complex64)
amps_3 = np.concatenate((amps_3, np.zeros(int(l/2+1)-len(amps_3))))

hLt3 = ht * np.conj(amps_3)
# hLt2 = ht_deg * np.conj(amps_2)
hL3 = make_transform_with_C_c2r(hLt3)
#%  load amps 4
l4 = 'NFW_2'
y4 = '002804'
M4 = '10-9'

af4 = '/home/paolo/Desktop/waveform/case_h/'+l4+'/amps/amps_case_h_H74_y'+y4+'_M'+M4+'_zL05'
amps_4 = np.fromfile(af4, np.complex64)
amps_4 = np.concatenate((amps_4, np.zeros(int(l/2+1)-len(amps_4))))

hLt4 = ht * np.conj(amps_4)
# hLt4 = ht_deg * np.conj(amps_4)
hL4 = make_transform_with_C_c2r(hLt4)
#%  load amps 5
l5 = 'NFW_2'
'''
y5 = '001061'

M5 = '10-9'

af5 = '/home/paolo/Desktop/waveform/case_h/'+l5+'/amps/amps_case_h_H74_y'+y5+'_M'+M5+'_zL05_norm'
amps_5 = np.fromfile(af5, np.complex64)
amps_5 = np.concatenate((amps_5, np.zeros(int(l/2+1)-len(amps_5))))

# hLt2 = ht * np.conj(amps_2)
'''
hLt5 = ht_deg * np.conj(amps_4)
hL5 = make_transform_with_C_c2r(hLt5)
#%% plot frequency domain

plt.figure()
plt.plot(fr, 2*fr*abs(ht), c='grey', label='unlensed')
plt.plot(fr, 2*fr*abs(ht_deg2), c='k', label='un - $z=0.8$')

plt.plot(fr, 2*fr*abs(hLt1),  c='red', label=l1)
plt.plot(fr, 2*fr*abs(hLt2), c='b', label=l2)
plt.plot(fr, 2*fr*abs(hLt3), c='green', label=lbl[2])
plt.plot(fr, 2*fr*abs(hLt4), dashes=(5, 5), c=df_clr[9], label=lbl[3]) #dashes=(5, 5)
plt.plot(fr, 2*fr*abs(hLt5), '-.', c='orangered', label=lbl[4])

# plt.plot(fr, 2*fr*abs(ht_deg), '--', c='k', linewidth=3, label='unlensed deg')

plt.xscale('log')
plt.yscale('log')
# plt.grid()
plt.xlim(9.9*10**-6, 1.1*10**-4)
plt.ylim(9.9*10**-17, 1.1*10**-15)
plt.xlabel('$f$ [Hz]', fontsize=fontsz)
plt.ylabel('characteristic strain $[\sqrt{4f^2|\\tilde{h}(f)|^2}]$', fontsize=fontsz)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)
# plt.legend(ncol=3,loc='lower center',bbox_to_anchor=(0.5,1.),fontsize=fontsz,framealpha=0.5)
plt.legend(ncol=2, loc='lower left',fontsize=fontsz-1,framealpha=0.1)
plt.show()
# %% save plot
plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/models_FD_ySIS1_M10-9.png',
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()





#%%
plt.figure()
[plt.hlines(i, 1,2, color=df_clr[n]) for n,i in enumerate(range(10))]
plt.show()




#%%

# PHASE!!

#%% calculate error of normalized phases
rho2 = 220
# rho1 = 10

d2 = 1/rho2

er1 = d2/ps_wfL1 + d2/ps_t
er2 = d2/ps_wfL2 + d2/ps_t
ea1 = (ps_wfL1/ps_t * er1)/2
ea2 = (ps_wfL2/ps_t * er2)/2
# %% plot phases with inset
iin = 0
# ifi = len(fr)#1262
norm = ps_t

# rho2 = 220
# rho1 = 10

d2 = 1/rho2/2
print(d2)
# d1 = 1/rho1/2
# print(d1)

c1 = 'blue'
c2 = 'orangered'

fig = plt.figure()
ax = fig.add_subplot(111)


# lensed 2
# '''
plt.fill_between(hpt.sample_frequencies, (ps_wfL2/norm)-ea2, (ps_wfL2/norm)+ea2, color=c2, alpha=0.5)#, label='$\\rho=%i$'%(rho2))
plt.plot(hpt.sample_frequencies, ps_wfL2/norm, c=c2, label=lbl[3])
# '''


# lensed 1

# '''
plt.fill_between(hpt.sample_frequencies, (ps_wfL1/norm)-ea1, (ps_wfL1/norm)+ea1, color=c1, alpha=0.3)#, label='$\\rho=%i$'%(rho2))
plt.plot(hpt.sample_frequencies, ps_wfL1/norm, c=c1, label=lbl[1])
# '''

plt.plot(hpt.sample_frequencies, ps_t/norm, linewidth=2, c='k', label='unlensed')
# plt.plot(hpt.sample_frequencies, ps_tdeg/norm, '--', linewidth=1,  c='r', label='un - $z=0.8$')

plt.xlim(1.3*10**-5, 1.1*10**-4)
# plt.ylim(0.979, 1.045)     # SIS
# plt.ylim(0.99845, 1.00045) # NFW-2
plt.ylim(0.9957, 1.0012)   #NFW

# plt.grid()
plt.xlabel('f [Hz]', fontsize = fontsz)
# plt.ylabel('$\phi/\phi_{unlensed}$')
plt.ylabel('$\phi$ [rad]', fontsize = fontsz)

plt.xscale('log')

plt.xticks(xt, labels=['$2\cdot10^{-5}$', '', '$4\cdot10^{-5}$','','$6\cdot10^{-5}$', '', '', '', '$10^{-4}$'])

plt.legend(ncol=4,loc='lower center',bbox_to_anchor=(0.5,1.),fontsize=fontsz,framealpha=0.5)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)


# inset
# ''' 
fontsz_i = fontsz-1
inset_axes_1 = inset_axes(ax, width="45%", height="45%", loc=4, 
                          bbox_to_anchor=(-0.03,0.02,1.,1.), bbox_transform=ax.transAxes,
                          axes_kwargs={'alpha':0.05, 'facecolor':'white'})
inset_axes_1.axes
plt.plot(hpt.sample_frequencies[1:], ps_1, c=c1, label=l1)
plt.plot(hpt.sample_frequencies[1:], ps_2, c=c2, label=l2)

plt.xlim(9.9*10**-6, 1.*10**-4)
# plt.ylim(-0.22, 0.23)   # SIS
# plt.ylim(-0.002, 0.008) # NFW-2
plt.ylim(-0.021, 0.023)#NFW

# plt.grid()
# plt.xlabel('f [Hz]', fontsize=fontsz_i)
plt.title('f [Hz]', fontsize=fontsz_i)
plt.ylabel('$\phi_{AF}$', fontsize = fontsz_i)
plt.xticks(xt, labels=['$2\cdot10^{-5}$', '', '$4\cdot10^{-5}$', '', '$6\cdot10^{-5}$', '', '', '', '$10^{-4}$'])
plt.xscale('log')
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz_i
                , bottom=False, labelbottom=False, top=True, labeltop=True)
# '''

plt.show()
# %% save plot
# plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/phase_shift_'+l1+'_we.png',
plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/phase_shift_'+l1+'vs'+l2+'.png',
            dpi=300, format='png', bbox_inches="tight")#, transparent=True
plt.close()
#%%

