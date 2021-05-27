#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 13 08:51:25 2021

@author: paolo
"""

#%%load packages and functions
from Dropbox.PhD.Python.program_py.packages.pyCBC_function import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
#%% make lens characteristic

case = 'd'

zS = 0.5
zL = 0.01
M_l = 500
dt = 10**-4#5*10**-5
y = 1.

DS = cosmo.luminosity_distance(zS).value #distance of the source in Mpc
MTOT = 100*(1+zS)
q = 1.
m1 = MTOT/(1+q)
m2 = MTOT-m1
hp, hc = get_td_waveform(approximant='IMRPhenomD', mass1=m1, mass2=m2,
                         delta_t=dt, f_lower=10, inclination=np.pi/3,
                         # s1z = 0.7, s2z = 0.2,
                         distance=DS, coa_phase=0.)

htp = make_transform_with_C_r2c(hp.data)
htc = make_transform_with_C_r2c(hc.data)
l = len(hp)
fr = np.arange(l/2+1)/(l*dt)
hp_t = htp*Fp(0.3, 0.4, 1.5)+htc*Fx(0.3, 0.4, 1.5)

amps_0 = np.array([AF_PM(i, M_l, zL, y, 1.) for i in fr[fr<700]])
amps_0 = np.concatenate((amps_0, np.zeros(len(hp_t)-len(amps_0))))
amps_0[0] = 0.
hp_L_t = hp_t * np.conj(amps_0)
#%% load hp
# uws = ['mt', 'q01', 's0702'] # case hez
uws = ['', '_q01T', '_s0702T'] # case d 
plt.figure()

for uw in uws:
    hp_t = np.fromfile('/home/paolo/Desktop/waveform/case_'+case+'/hp_fd_H74'+uw)
    l = len(hp_t)-2
    hp_t = unite(hp_t)
    fr = np.arange(l/2+1)/(l*dt)
    plt.plot(fr, 2*fr*abs(hp_t), label=uw[1:-1])

# plt.xlim(29, 816) # case hez
# plt.ylim(10**-20, 1.5*10**-17) # case hez
plt.xlim(9.5, 460) # case h
plt.ylim(2.*10**-19, 7.*10**-17) # case d
plt.grid()
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.xlabel('f [Hz]')
plt.ylabel('characteristic strain $[\sqrt{4f^2|\\tilde{h}(f)|^2}]$')
'''
pylab.figure()
pylab.plot(fr, 2*fr*abs(hp_t), c='grey', label='unlensed')

pylab.xscale('log')
pylab.yscale('log')
pylab.grid()
pylab.legend()
# '''
#%%
directory = '/home/paolo/Desktop/waveform/case_'+case+'/'+lens+'/amps/'
# uws = ['mt', 'q01', 's0702'] # case hez
uws = ['', '_q01T', '_s0702T'] # case d
uw = uws[0]
file_td = directory[:-8]+'hp_fd_H74'+uw
# f_a = directory+'amps_case_hez_H74_y2_M10-4_zL01' #_q01
mass = '100'
f_a = directory+'amps_case_d_H74_y1_M'+mass+'_zL001' #_q01

hp_L_t = make_hpL(f_a, file_td, fd = True)
print(uw)

hp_t = np.fromfile(file_td)
l = len(hp_t)-2
hp_t = unite(hp_t)
fr = np.arange(l/2+1)/(l*dt)
#%%
lamdas = ['05', '075', '095', '11']
clrs = ['darkblue', 'cyan', 'green', 'orange']

r = 0.52

fig = plt.figure()
ax = fig.add_subplot(111)

for n, la in enumerate(lamdas):
    lMSD = float(la[:1]+'.'+la[1:])
    print(lMSD)
    amps_1 = np.array([AF_PM(i, M_l, zL, y, lMSD) for i in fr[fr<700]])
    amps_1 = np.concatenate((amps_1, np.zeros(len(hp_t)-len(amps_1))))
    amps_1[0] = 0.
    hp_t_1 = hp_t * np.conj(amps_1)
    plt.plot(fr, 2*fr*abs(hp_t_1), c =clrs[n], label='$\lambda = $' + la[:1] + '.' + la[1:])
    '''
    file_amps = f_a + '_L' + la
    hp_L_t_l = make_hpL(file_amps, file_td, fd = True)
    plt.plot(fr, 2*fr*abs(hp_L_t_l), c =clrs[n], label='$\lambda = $' + la[:1] + '.' + la[1:])
    
    if la == '05':
        hp_L_05 = hp_L_t_l
    '''
    
plt.plot(fr, 2*fr*abs(hp_L_t), c = 'red', label='$\lambda = 1$')
plt.plot(fr, 2*fr*abs(hp_t), c = 'grey', label='unlensed')

# plt.xlim(9.8, 460)
# plt.ylim(3.*10**-19, 2.*10**-16)
plt.xlim(9.8, 200)
plt.ylim(3.*10**-19, 2.*10**-17)
plt.grid()

# plt.legend(fontsize=8)
handles, labels = plt.gca().get_legend_handles_labels()
# order = [5, 0, 1, 2, 4, 6, 3]
# order = [5, 2, 7, 0, 4, 3, 1, 6]
order = [5, 2, 0, 4, 1, 3]
plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=fontsz, loc='lower left')

plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)
# plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz, labelright=True, right=True, labelleft=False, left=False)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$f$ [Hz]', fontsize=fontsz)
plt.ylabel('characteristic strain $[\sqrt{4f^2|\\tilde{h}(f)|^2}]$', fontsize=fontsz)

# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], 
           # ncol=3,loc='lower center',bbox_to_anchor=(0.5,1.),fontsize=fontsz,framealpha=0.5)

'''
inset_axes_1 = inset_axes(ax, width="55%", height="55%", loc=3, axes_kwargs={'alpha':0.05, 'facecolor':'white'})
# inset_axes_1.axes
plt.plot(fr, 2*fr*abs(hp_L_t), c = 'red')#, label='$\lambda = 1$')
plt.plot(fr, 2*fr*abs(hp_L_05), c =clrs[0])#, label='$\lambda = 0.5$')
plt.plot(fr, 2*fr*abs(hp_L_t)/r, c = 'k', label='$\lambda = 1$ - rescaled')
plt.xscale('log')
plt.yscale('log')
plt.xlim(9.8, 300)
plt.ylim(4.*10**-17, 1.5*10**-16)
# plt.ylim(4.*10**-17, 0.7*10**-16)
plt.grid()
plt.legend(loc='lower center',bbox_to_anchor=(0.38,.97),fontsize=fontsz,framealpha=0.5)
plt.tick_params(axis='both',which='both',left=False, bottom=False, labelbottom=False, labelleft=False)
# '''

plt.show()
#%%
# plt.savefig('/home/paolo/Desktop/waveform/MSD/plots/case_'+case+'_M'+mass+'_y1_fd'+uw+'_sp.png', 
plt.savefig('/home/paolo/Dropbox/PhD/plots/FD/q1_TD_M500_y1_zs05.png',
            dpi=300, format='png', bbox_inches="tight") #, transparent=True)
plt.close()
#%%plot time domain
hp = np.fromfile('/home/paolo/Desktop/waveform/case_d/hp_td_H74_s0702')
hp_L = make_hpL('/home/paolo/Desktop/waveform/case_d/PM/amps/amps_case_d_H74_y5_M100_zL001',
                '/home/paolo/Desktop/waveform/case_d/hp_fd_H74_s0702T') #q01T s0702T
t = np.arange(0, len(hp)*dt, dt)
tL = np.arange(0, len(hp_L)*dt, dt)
#%%plot
pu = np.where(hp==min(hp))[0][0]
pl = np.where(hp_L==min(hp_L))[0][0]

ii = 0
ff = len(hp)

plt.figure()
# plt.plot(t[ii:ff], hp[ii:ff], c='grey', label='unlensed')
# plt.plot(t, np.roll(hp, 78630-4189820), c='grey', label='unlensed')
plt.plot(t, np.roll(hp, pl-pu), c='grey', label='unlensed')

plt.plot(tL, hp_L, c='orangered', label='lensed')
# plt.plot(tL, np.roll(hp_L, pu-pl), c='orangered', label='lensed')

plt.grid()
plt.legend(loc='lower left')
# plt.xlim(3.75, 4.025)
plt.xlim(3.75, 4.05)
plt.xlabel('t [s]')
plt.ylabel('strain')
plt.show()
#%%
plt.savefig('/home/paolo/Desktop/waveform/MSD/plots/len_unlens_y5_M100_TAB.png', dpi=300, format='png', bbox_inches="tight")
#%% inner product

#%% S/N calculation
# 2008.12814
#  load files
dt = 5*10**-5
# Sn = pd.read_csv('/home/paolo/Desktop/aligo_O3actual_L1.txt', sep=' ')
# Sn = pd.read_csv('/home/paolo/Desktop/AplusDesign.txt', sep='  ', header=None)
Sn = pd.read_csv('/home/paolo/Downloads/curves_Sep_2018/aplus.txt', sep=' ', header=None)
Sn.columns = ['fr', 'cs']
# hT = np.fromfile('/home/paolo/Desktop/waveform/case_d/hp_fd_H74')
# file_t = '/home/paolo/Desktop/waveform/case_d/hp_fd_H74_s0702T'
file_t = '/home/paolo/Desktop/waveform/case_d/hp_fd_H74_s0702T'
file_a = '/home/paolo/Desktop/waveform/case_d/PM/amps/amps_case_d_H74_y05_M500_zL001'

hT = make_hpL(file_a, file_t, fd=True)

l = len(hT)*2
# hT = unite(hT)
fr = np.arange(l/2)/(l*dt)
lam = '11'
# hL = make_hpL('/home/paolo/Desktop/waveform/case_d/PM/amps/amps_case_d_H74_y1_M500_zL001_L'+lam,
              # '/home/paolo/Desktop/waveform/case_d/hp_fd_H74', fd=True)
#%%
'''
#  https://stackoverflow.com/questions/44811581/scipy-how-to-integrate-a-linearly-interpolated-function
# CHECK! Probabilmente c'è qualcosa che non va
snr_1 = inner_p(hL, fr, hT, fr, Sn.cs**2, Sn.fr)
snr_1_T = inner_p(hT, fr, hT, fr, Sn.cs**2, Sn.fr)
snr_2 = inner_p(hT, fr, hT, fr, Sn.cs**2, Sn.fr)
snr = snr_1/np.sqrt(snr_2)
snr_T = snr_1_T/np.sqrt(snr_2)
print('S/N%s = %.3e\nS/N template = %.3e\nratio%s= %.4f'%(' '*9, snr, snr_T, ' '*8, snr/snr_T))
'''

#%%
case = 'd'
ys = '1'
ms = '100'
dt = 5*10**-5
f_min = 10
f_max = 500#175#155

lsa = [[],[],[]]
snrs = [[],[],[]]

# Sn = pd.read_csv('/home/paolo/Desktop/aligo_O3actual_L1.txt', sep=' ',  header=None)
# Sn = pd.read_csv('/home/paolo/Desktop/AplusDesign.txt', sep='  ',  header=None)
Sn = pd.read_csv('/home/paolo/Downloads/curves_Sep_2018/aplus.txt', sep=' ', header=None)
Sn.columns = ['fr', 'cs']
for n, Temp in  enumerate(['', '_q01T', '_s0702T']):
# for n, Temp in  enumerate(['', '_q01', '_s0702']):
    print(n+1)
    file_t = '/home/paolo/Desktop/waveform/case_'+case+'/hp_fd_H74'+Temp
    file_a = '/home/paolo/Desktop/waveform/case_'+case+'/PM/amps/amps_case_'+case+'_H74_y'+ys+'_M'+ms+'_zL001'
    hT = make_hpL(file_a, file_t, fd=True)
    l = len(hT)*2
    fr = np.arange(l/2)/(l*dt)
    snr_1_T = inner_p(hT, fr, hT, fr, Sn.cs**2, Sn.fr, f_min, f_max)
    snr_T = np.sqrt(snr_1_T)

    for lam in ['04', '05', '06', '07', '075', '08', '09', '095', '1', '105', '11', '12']: 
    # for lam in ['05', '075', '095', '11'] :
        hL = make_hpL(file_a+'_L'+lam, file_t, fd=True)
        snr_1 = inner_p(hT, fr, hL, fr, Sn.cs**2, Sn.fr, f_min, f_max)
        snr_d = inner_p(hL, fr, hL, fr, Sn.cs**2, Sn.fr, f_min, f_max)
        snr = snr_1/np.sqrt(snr_d)
        llam = llam = lam[:1] + '.' + lam[1:]
        snrs[n].append(snr/snr_T)
        lsa[n].append(float(llam))
#%
ttls = ['$q=1$', '$q=0.1$', '$s_{1,2;z} = [0.7, 0.2]$']
colrs = ['forestgreen', 'b', 'red']
plt.figure()
for i in range(3):
    plt.plot(lsa[i], snrs[i], linewidth=2, c=colrs[i], label=ttls[i]) 
    # plt.scatter(lsa[i], snrs[i], c=colrs[i], marker='.') 

plt.scatter(1, 1, label='template', c='k', marker='s')
# plt.xlim(0.35, 1.15)
plt.ylim(0.5, 1.01)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)

plt.xlabel('$\lambda$', fontsize=fontsz)
plt.ylabel(r'$\rho/\rho_{opt}$', fontsize=fontsz)
plt.grid()
# plt.legend()
plt.show()
#%%
plt.savefig('/home/paolo/Dropbox/PhD/plots/SNR_M'+ms+'_y'+ys+'.png', dpi=300, format='png', transparent=True, bbox_inches="tight")

# plt.savefig('/home/paolo/Desktop/waveform/MSD/plots/SNR_M500_y05_T.png', dpi=300, format='png', bbox_inches="tight")

#%%
plt.figure()
plt.plot(fr, 2*fr*abs(hT), label='template - S/N = %.2e'%(snr_T))
llam = lam[:1] + '.' + lam[1:]
plt.plot(fr, 2*fr*abs(hL), label='MST $\lambda = %s$ - S/N = %.2e'%(llam, snr))
plt.plot(Sn.fr, Sn.cs**2, c='k')
# plt.xlim(9, 200)
# plt.ylim(5*10**-18, 2*10**-16)
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
#%%

#%%
now = datetime.now()

case = 'd'
yy = '1'
ms = '500'
dt = 5*10**-5
f_min = 10
f_max = 500#155#175#155

# ys = ['005', '05', '07', '1', '2']
ys = np.arange(0., 5.1, 0.1)
print(len(ys))
lsa = [[],[],[]]
snrs = [[],[],[]]

# Sn = pd.read_csv('/home/paolo/Desktop/aligo_O3actual_L1.txt', sep=' ',  header=None)
# Sn = pd.read_csv('/home/paolo/Desktop/AplusDesign.txt', sep='  ',  header=None)
Sn = pd.read_csv('/home/paolo/Downloads/curves_Sep_2018/aplus.txt', sep=' ', header=None)
Sn.columns = ['fr', 'cs']
for n, Temp in  enumerate(['', '_q01T', '_s0702T']):
# for n, Temp in  enumerate(['', '_q01', '_s0702']):
    print(n+1)
    file_t = '/home/paolo/Desktop/waveform/case_'+case+'/hp_fd_H74'+Temp
    file_a = '/home/paolo/Desktop/waveform/case_'+case+'/PM/amps/amps_case_'+case+'_H74_y'+yy+'_M'+ms+'_zL001'
    hT = make_hpL(file_a, file_t, fd=True)
    l = len(hT)*2
    fr = np.arange(l/2)/(l*dt)
    snr_1_T = inner_p(hT, fr, hT, fr, Sn.cs**2, Sn.fr, f_min, f_max)
    snr_T = np.sqrt(snr_1_T)

    for y in ys: 
    # for lam in ['05', '075', '095', '11'] :
        '''
        file_a1 = '/home/paolo/Desktop/waveform/case_'+case+'/PM/amps/amps_case_'+case+'_H74_y'+y+'_M'+ms+'_zL001'
        hL = make_hpL(file_a1, file_t, fd=True)
        '''
        # yf = float(y[:1]+'.'+y[1:])
        amps = np.array([AF_PM(i, float(ms), zL, y, 1.) for i in fr[fr<700]])
        amps[0] = 0.
        hL = make_hpL(amps, file_t, fd=True)
        # '''
        snr_1 = inner_p(hT, fr, hL, fr, Sn.cs**2, Sn.fr, f_min, f_max)
        snr_d = inner_p(hL, fr, hL, fr, Sn.cs**2, Sn.fr, f_min, f_max)
        snr = snr_1/np.sqrt(snr_d)
        llam = llam = lam[:1] + '.' + lam[1:]
        snrs[n].append(snr/snr_T)
        lsa[n].append(y)

print('\ntime: \n', datetime.now()-now)    
#%%
ttls = ['$q=1$', '$q=0.1$', '$s_{1,2;z} = [0.7, 0.2]$']
colrs = ['forestgreen', 'b', 'red']
plt.figure()
for i in range(3):
    plt.plot(lsa[i], snrs[i], linewidth=2, c=colrs[i], label=ttls[i]) 
    # plt.scatter(lsa[i], snrs[i], c=colrs[i], marker='.') 

plt.scatter(1, 1, label='template', c='k', marker='s')
# plt.xlim(0.35, 1.15)
# plt.ylim(0.5, 1.01)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)

plt.xlabel('$y$', fontsize=fontsz)
plt.ylabel(r'$\rho/\rho_{opt}$', fontsize=fontsz)
plt.grid()
plt.legend()
plt.show()

#%%
plt.savefig('/home/paolo/Desktop/waveform/new/SNR_M'+ms+'_ys_TF500.png', dpi=300, format='png', transparent=True, bbox_inches="tight")
#%%
amps1 = np.array([AF_PM(i, float(ms), zL, 0.2, 1.) for i in fr[fr<700]])
amps1[0] = 0.
hL1 = make_hpL(amps1, '/home/paolo/Desktop/waveform/case_'+case+'/hp_fd_H74')#, fd=True)

amps2 = np.array([AF_PM(i, float(ms), zL, 2.4, 1.) for i in fr[fr<700]])
amps2[0] = 0.
hL2 = make_hpL(amps2, '/home/paolo/Desktop/waveform/case_'+case+'/hp_fd_H74')#, fd=True)

#%%
plt.figure()
'''
plt.plot(fr, 2*fr*abs(hL1), label='y=0.1')
plt.plot(fr, 2*fr*abs(hL2), label='y=2.7')
plt.xscale('log')
plt.yscale('log')
# '''
plt.plot(hL1, label='y=0.2')
plt.plot(hL2, label='y=2.4')
# '''
plt.legend()
plt.show()


