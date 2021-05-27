#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:59:35 2021

@author: paolo
"""

#%%load packages and functions
from Dropbox.PhD.Python.program_py.packages.pyCBC_function import *
from Dropbox.PhD.Python.program_py.packages import LISA as li
lisa = li.LISA()
#%% define case
case = 'h'#'d'
df = 10**-6
dt = 1
zS = 1.
DS = cosmo.luminosity_distance(zS).value #distance of the source in Mpc
zL = 0.5
Dl = D_l(zL, D_a=True)
DL = cosmo.angular_diameter_distance(zL).value*pctomt*10**6
Ml = 10**9
f_min = 10**-5
MTOT = 10**8*(1+zS)
q = 1.
m1 = MTOT/(1+q)
m2 = MTOT-m1
#%% frequency domain
'''
ht, htc = get_fd_waveform(approximant='IMRPhenomPv3', mass1=m1, mass2=m2,
                         delta_f=df, f_lower=f_min, inclination=np.pi/3,
                         # s1z = 0.7, s2z = 0.2,
                         distance=DS, coa_phase=np.pi/4)

fr = ht.sample_frequencies
# ht = hp.data*Fp(0.3, 0.4, 1.5)+hc.data*Fx(0.3, 0.4, 1.5)
'''
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
#%% degenerate wf
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
#'''
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

#%%  load amps 2
l2 = 'NFW_2'
# y2 = '001061'
y2 = '002804'
M2 = '10-9'

af2 = '/home/paolo/Desktop/waveform/case_h/'+l2+'/amps/amps_case_h_H74_y'+y2+'_M'+M2+'_zL05_norm'
amps_2 = np.fromfile(af2, np.complex64)
amps_2 = np.concatenate((amps_2, np.zeros(int(l/2+1)-len(amps_2))))

# hLt2 = ht * np.conj(amps_2)
hLt2 = ht_deg * np.conj(amps_2)
hL2 = make_transform_with_C_c2r(hLt2)
#%%  load amps 3
l3 = 'GNFW_2'
y3 = '00634'
M3 = '10-9'

af3 = '/home/paolo/Desktop/waveform/case_h/'+l3+'/amps/amps_case_h_H74_y'+y3+'_M'+M3+'_zL05'#_FP'
amps_3 = np.fromfile(af3, np.complex64)
amps_3 = np.concatenate((amps_3, np.zeros(int(l/2+1)-len(amps_3))))

hLt3 = ht * np.conj(amps_3)
# hLt2 = ht_deg * np.conj(amps_2)
hL3 = make_transform_with_C_c2r(hLt3)
#%% plot frequency domain
plt.figure()
plt.plot(fr, 2*fr*abs(ht), c='grey', label='unlensed')
plt.plot(fr, 2*fr*abs(hLt1),  c='red', label=l1)#+' - y = '+y1[:1]+'.'+y1[1:])#+' - zL = '+zL1[:1]+'.'+zL1[1:])
# plt.plot(fr, 2*fr*abs(hLt2), c='orangered', label=l2+ ' - $z_s=0.55$')
plt.plot(fr, 2*fr*abs(hLt2), c='b', label=l2)#+' - y = '+y2[:1]+'.'+y2[1:]+' - M='+M2)
plt.plot(fr, 2*fr*abs(hLt3), c='green', label=l3)#+' - y = '+y3[:1]+'.'+y3[1:])#+' - M='+M3)

# plt.plot(fr, 2*fr*abs(ht_deg), '--', c='k', linewidth=3, label='unlensed deg')

plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.xlim(9.9*10**-6, 1.1*10**-4)
plt.ylim(9.9*10**-17, 1.1*10**-15)
plt.xlabel('$f$ [Hz]', fontsize=fontsz)
plt.ylabel('characteristic strain $[\sqrt{4f^2|\\tilde{h}(f)|^2}]$', fontsize=fontsz)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)
plt.legend(ncol=3,loc='lower center',bbox_to_anchor=(0.5,1.),fontsize=fontsz,framealpha=0.5)
plt.show()
# %% save plot
plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/models_FD_ySIS1_M10-9.png',
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()

#%% find peaks 
'''
puv = find_peaks(hp)[0]
pu = max(hp[puv])
pu_ind = np.where(hp==pu)[0][0]
'''
pu = np.where(hp.data==max(hp.data))[0][0]
p1v = np.where(hL1==max(hL1))[0][0]
p2v = np.where(hL2==max(hL2))[0][0]
#%% plot time domain
tsc = 10**5
ssc = 10**-16
plt.figure()

# plt.plot(hp.sample_times/tsc, hp/ssc, c='grey', label='unlensed')
plt.plot(hp.sample_times/tsc, np.roll(hL1/ssc, pu-p1v), c='red', label='SIS')
plt.plot(hp.sample_times/tsc, np.roll(hL2/ssc, pu-p2v), c='b', label=l2)
# plt.plot(hp.sample_times/tsc, hL2/ssc, c='b', label=l2)
# plt.plot(hp_deg.sample_times/tsc, hp_deg/ssc, c='k', linewidth=1, label='unlensed deg - coa:$2.2\pi/4$')
             
plt.grid()
# plt.xlim(-2.5, 0.8)
# plt.ylim()
plt.xlabel('$t$ [$10^5$ s]', fontsize=fontsz)
plt.ylabel('strain [$10^{-16}$]', fontsize=fontsz)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)
plt.legend(ncol=3,loc='lower center',bbox_to_anchor=(0.5,1.),fontsize=fontsz,framealpha=0.5)
plt.show()
# %% save plot
plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/SISvsNFW_y1_M10-9_TD.png',
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()




#%% make SNR
f_max = vf_fin(f_min, m1, m2)
'''
dd = 'o3_l1'#'o3_l1'#'aligo'#'aplus'
Sn = pd.read_csv('/home/paolo/Downloads/curves_Jan_2020/'+dd+'.txt', sep=' ', header=None)
Sn.columns = ['fr', 'cs']
'''

fd  = np.logspace(np.log10(1.0e-5), np.log10(1.0e0), 1000)
Sn = lisa.Sn(fd)

hS = hLt2.data
hT = hLt1.data
# hT = ht.data
# hT = ht_deg.copy()
#%%
pwr = 1
snr_n = inner_p(hS, fr, hT, fr, Sn**pwr, fd, f_min, f_max)
snr_d = inner_p(hT, fr, hT, fr, Sn**pwr, fd, f_min, f_max)
SNR_T = np.sqrt(np.real(inner_p(hS, fr, hS, fr, Sn**pwr, fd, f_min, f_max)))

sg = np.sqrt(11.8)#3
th = 1-(1/2*(sg/SNR_T)**2)


print('detected SNR = ',np.real((snr_n/np.sqrt(snr_d))))
print('optimal SNR  = ',np.real(SNR_T))
print('template SNR = ',np.real(np.sqrt(snr_d)))
print('ratio SNR    = ',np.real((snr_n/np.sqrt(snr_d)))/np.real(SNR_T))
#%%
plt.figure()

plt.plot(fd, Sn**1, label='LISA')
plt.plot(fr, 2*fr*abs(ht), c='grey', label='unlensed')
plt.plot(fr, 2*fr*abs(hLt1), c='red', label='SIS')
plt.plot(fr, 2*fr*abs(hLt2), c='b', label=l2)

plt.xscale('log')
plt.yscale('log')
plt.grid()
plt.legend()

plt.show()


