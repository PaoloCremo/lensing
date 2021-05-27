#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 30 11:25:30 2021

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
Ml = 10**9
f_min = 10**-5
MTOT = 10**8*(1+zS)
q = 1.
m1 = MTOT/(1+q)
m2 = MTOT-m1
#%% frequency domain
# '''
hpt, hct = get_fd_waveform(approximant='IMRPhenomPv3', mass1=m1, mass2=m2,
                         delta_f=df, f_lower=f_min, inclination=np.pi/3,
                         # s1z = 0.7, s2z = 0.2,
                         distance=DS, coa_phase=0)#np.pi/2)

fr = hpt.sample_frequencies
#%% deg frequency domain
zS_deg = 0.5
DS_deg = cosmo.luminosity_distance(zS_deg).value 
hpt_deg, hct_deg = get_fd_waveform(approximant='IMRPhenomPv3', mass1=m1, mass2=m2,
                         delta_f=df, f_lower=f_min, inclination=np.pi/3,
                         # s1z = 0.7, s2z = 0.2,
                         distance=DS_deg, coa_phase=0)

#%% load amps FP
l1 = 'NFW'
y1 = '001061'
l2 = 'NFW_2'
y2 = '002804'

afFP1 = '/home/paolo/Desktop/waveform/case_h/'+l1+'/amps/amps_case_h_H74_y'+y1+'_M10-9_zL05_FP'
ampsFP_1 = np.fromfile(afFP1, np.complex64)
hLtFP_1 = hpt * np.conj(ampsFP_1)

afFP2 = '/home/paolo/Desktop/waveform/case_h/'+l2+'/amps/amps_case_h_H74_y'+y2+'_M10-9_zS05_zL025_FP' #y002804
ampsFP_2 = np.fromfile(afFP2, np.complex64)
# hLtFP_2 = hpt * np.conj(ampsFP_2)
hLtFP_2 = hpt_deg * np.conj(ampsFP_2)

#%% plot frequency domain
plt.figure()
# plt.plot(fr, 2*fr*abs(ht_deg), c='k', linewidth=1, label='un - $z=0.8$')
plt.plot(fr, 2*fr*abs(hpt), c='grey', label='un - $z=1$')
plt.plot(fr, 2*fr*abs(hLtFP_1),  c='red', label=l1)
# plt.plot(hpt.sample_frequencies, 2*hpt.sample_frequencies*abs(hLtFP_1))
plt.plot(fr, 2*fr*abs(hLtFP_2),  c='b', label=l2)#+' - $z_S=%.2f$'%(zS_deg)) #dashes=(10, 10),
# plt.plot(hpt.sample_frequencies, 2*hpt.sample_frequencies*abs(hLtFP_2))

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
# plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/SISNFW_y1_M10-9_TD.png',
plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/NFWvsNFW2_pshift.png',
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()

#%% phase shift study
af1 = ampsFP_1
af2 = ampsFP_2
ps_1 = []
for i in range(1, len(af1)): #hpt
    ps_1.append(- np.real(1j*np.log(af1[i]/np.abs(af1[i]))))

ps_2 = []
for i in range(1, len(af2)):
    ps_2.append(- np.real(1j*np.log(af2[i]/np.abs(af2[i]))))
   
ps_1 = np.array(ps_1)
ps_2 = np.array(ps_2)

xt = [a*10**-5 for a in range(2, 11)]
# xt = [a*10**-5 for a in range(2, 11, 2)]
#%% plot phase AF

rhoAF = 220
# rho1 = 10

dAF = 1/rhoAF/2
print(round(dAF, 5))

plt.figure()
plt.plot(hpt.sample_frequencies[1:], ps_1, c='red', label=l1) #np.unwrap(ps_1)
plt.plot(hpt.sample_frequencies[1:], np.array(ps_2), c=df_clr[9], label=l2)

# plt.fill_between(hpt.sample_frequencies[1:], ps_1-dAF, ps_1+dAF, color='paleturquoise', alpha=0.5)
# plt.fill_between(hpt.sample_frequencies[1:], ps_2-dAF, ps_2+dAF, color='orangered', alpha=0.5)

plt.xlim(9.9*10**-6, 1.1*10**-4)
# plt.ylim(-0.2, 30.5)
plt.grid()
plt.xlabel('f [Hz]', fontsize=fontsz)
# plt.ylabel('phase shift [$\pi$]', fontsize = fontsz)
plt.ylabel('$\phi_{AF}$', fontsize = fontsz)
plt.xticks(xt, labels=['$2\cdot10^{-5}$', '', '$4\cdot10^{-5}$', '', '$6\cdot10^{-5}$', '', '', '', '$10^{-4}$'])
plt.xscale('log')
# plt.legend(ncol=2,loc='lower center',bbox_to_anchor=(0.5,1.),fontsize=fontsz,framealpha=0.5)
plt.legend(fontsize=fontsz, framealpha=0.5)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)
plt.show()
# %% save plot
plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/phase_shift_'+l2+'.png',
# plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/'+l1+'vs'+l2+'_pshift.png',
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()

#%%calculate phases

# ps_wf = [cmath.phase(x) for x in ht]
ps_wfL1 = phase_from_frequencyseries(hLtFP_1)#*10**10)
ps_wfL2 = phase_from_frequencyseries(hLtFP_2)#*10**10)
ps_t = phase_from_frequencyseries(hpt)
# ps_tdeg = phase_from_frequencyseries(hpt_deg)
# ps_t2 = phase_from_frequencyseries(hpt2)

#%% plot phases
iin = 0
# ifi = len(fr)#1262
norm = 1#ps_t

rho2 = 100
# rho1 = 10

d2 = 1/rho2/2
print(d2)
# d1 = 1/rho1/2
# print(d1)

plt.figure()


# lensed 2
# '''
plt.fill_between(hpt.sample_frequencies, (ps_wfL2-d2)/norm, (ps_wfL2+d2)/norm, color='paleturquoise', alpha=0.5, label='$\\rho=%i$'%(rho2))
plt.plot(hpt.sample_frequencies, ps_wfL2/norm, c=df_clr[9], label=l2)
# '''

# lensed 1
# '''
plt.fill_between(hpt.sample_frequencies, (ps_wfL1-d2)/norm, (ps_wfL1+d2)/norm, color='orangered', alpha=0.5, label='$\\rho=%i$'%(rho2))
plt.plot(hpt.sample_frequencies, ps_wfL1/norm, c='r', label=l1)
# '''

plt.plot(hpt.sample_frequencies, ps_t/norm, linewidth=1, c='k', label='unlensed')
# plt.plot(hpt.sample_frequencies, ps_tdeg/norm, '--', linewidth=1,  c='r', label='un - $z=0.8$')

plt.xlim(9.9*10**-6, 3*10**-4)

# plt.xlim(1.75*10**-5, 1.05*10**-4)
# plt.ylim(0.9992, 1.0001)

plt.grid()
plt.xlabel('f [Hz]', fontsize = fontsz)
# plt.ylabel('$\phi/\phi_{unlensed}$')
plt.ylabel('$\phi$ [rad]', fontsize = fontsz)

plt.xscale('log')

plt.legend(ncol=4,loc='lower center',bbox_to_anchor=(0.5,1.),fontsize=fontsz,framealpha=0.5)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)
plt.show()
# %% save plot
# plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/shifted_phase_SISvs'+l2+'.png',
# plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/SIS_shifted_phase_norm.png', 
plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/'+l2+'_shifted_phase_norm.png',            
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()
#%% calculate error of normalized phases
rho2 = 100
# rho1 = 10

d2 = 1/rho2

er1 = d2/ps_wfL1 + d2/ps_t
er2 = d2/ps_wfL2 + d2/ps_t
ea1 = (ps_wfL1/ps_t * er1)/2
ea2 = (ps_wfL2/ps_t * er2)/2


#%% plot residuals
norm = ps_t

print(d2)
# d1 = 1/rho1/2
# print(d1)

plt.figure()

#fill - paleturquoise
plt.fill_between(hpt.sample_frequencies, (ps_wfL2/norm)-ea2, (ps_wfL2/norm)+ea2, color='green', alpha=0.3)#, label='$\\rho=%i$'%(rho2))
plt.plot(hpt.sample_frequencies, ps_wfL2/norm, c='green', label=l2)
# lensed 2
plt.fill_between(hpt.sample_frequencies, (ps_wfL1/norm)-ea1, (ps_wfL1/norm)+ea1, color='red', alpha=0.3)#, label='$\\rho=%i$'%(rho2))
plt.plot(hpt.sample_frequencies, ps_wfL1/norm, c='red', label=l1)

# unlensed
# plt.plot(hpt.sample_frequencies, ps_t/norm, linewidth=1, c='k', label='un - $z=1$')
# plt.plot(hpt.sample_frequencies, ps_tdeg/norm, '--', linewidth=1,  c='r', label='un - $z=0.8$')

plt.xlim(1.3*10**-5, 1*10**-4)
# plt.ylim(0.9958, 1.0012)

plt.grid()
plt.xlabel('f [Hz]', fontsize = fontsz)
# plt.ylabel('$\phi/\phi_{unlensed}$')
plt.ylabel('phase ratio', fontsize = fontsz)

plt.xscale('log')
plt.xticks([], labels=[])
plt.xticks([2*10**-5, 3*10**-5, 4*10**-5, 6*10**-5, 10**-4], 
           labels=['$2\cdot10^{-5}$', '', '$4\cdot10^{-5}$', '$6\cdot10^{-5}$', '$10^{-4}$'])
# plt.legend(ncol=4,loc='lower center',bbox_to_anchor=(0.5,1.),fontsize=fontsz,framealpha=0.5)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)
plt.show()
# %% save plot
plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/shifted_phase_'+l1+'vs'+l2+'_norm.png',
# plt.savefig('/home/paolo/Dropbox/PhD/plots/lens_masses/'+l2+'_shifted_phase_norm.png',            
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()
#%%

def resize_amp(amp, fr, fr_new):
    rp = interp1d(fr, np.real(amp))
    ip = interp1d(fr, np.imag(amp))
    # size = 2*size
    # df = fr[1]
    # u = df*size/2
    # dt = size+1/(size)*1/(2*u)
    # fr_new = np.arange(size/2)/(size*dt)
    
    # amp_new = [np.complex(rp(x), ip(x)) for x in fr_new] 
    amp_new = []
    i = 1/20
    for n,x in enumerate(fr_new):
        amp_new.append(np.complex(rp(x), ip(x)))
        if n>i*len(fr_new):
            print('%i%% done'%(i*100))
            i += 1/20
        
        
    return amp_new

def transf_f (f):
    ome = 8*np.pi*G*Ml*(1+zL)*smtokg/c**3*f
    
    return ome




















