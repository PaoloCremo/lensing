#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 26 11:20:58 2021

@author: paolo
"""

#%%load packages and functions
# from Dropbox.PhD.Python.program_py.packages.pyCBC_function import *
from packages.pyCBC_function import *
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
dc = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
              '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
              '#bcbd22', '#17becf']
#%%
zS_WO = zS_INT = 0.1
zS_GO = 0.5

zL_WO = zL_INT = 0.01
zL_GO = 0.1

Ml_INT = 500
Ml_WO = 100
Ml_GO = 10**4

y_INT = y_WO = 1.
y_GO = 5.

DS = cosmo.luminosity_distance(zS_INT).value #distance of the source in Mpc


MTOT_WO = MTOT_INT = 100*(1+zS_INT)
MTOT_GO = 60*(1+zS_GO)

q = 1.
# m1 = MTOT/(1+q)
# m2 = MTOT-m1
#%%
df = 0.2
hp_INT, hc_INT = get_fd_waveform(approximant='IMRPhenomHM', mass1=MTOT_INT/(1+q), 
                                 mass2=MTOT_INT*q/(q+1),
                                 delta_f=df, f_lower=10, inclination=np.pi/3,
                                 distance=DS, coa_phase=0.)

fr = hp_INT.sample_frequencies.data
hp_t_INT = hp_INT.data*Fp(0.3, 0.4, 1.5)+hc_INT.data*Fx(0.3, 0.4, 1.5)

#%%
Ml = Ml_INT
zL = zL_INT
y = y_INT
lambdas = np.arange(0.4, 1.21, 0.01)
lambdas_i = np.arange(0., 1.21, 0.005)
print(len(lambdas_i))
dd = 'o3_l1'#'o3_l1'#'aligo'#'aplus'
Sn = pd.read_csv('/home/paolo/Downloads/curves_Jan_2020/'+dd+'.txt', sep=' ', header=None)
Sn.columns = ['fr', 'cs']
pwr = 2


#%% INT
amps = pd.DataFrame(columns=lambdas_i.round(3))
for n,lam in enumerate(lambdas_i.round(3)):
    a = []
    for f in fr:
        a.append(AF_PM(f, Ml, zL, y, lam))
    a[0] = 0.
    amps[lam] = a
    print('\r%.1f %% done'%(n/len(lambdas_i)*100), end='')

f_min = 11
f_max = vf_fin(f_min, MTOT_INT/(1+q), MTOT_INT*q/(q+1))

print('\nf max = %.2f'%(f_max))

snrs = []
hS = hp_t_INT * np.conj(amps[1.])
for lam in lambdas_i.round(3):
    hT = hp_t_INT * np.conj(amps[lam])#make_hpL(amps[lam], hp_t_INT, fd=True)
    snr_n = inner_p(hS, fr, hT, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)
    snr_d = inner_p(hT, fr, hT, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)
    SNR = np.real((snr_n/np.sqrt(snr_d)))
    snrs.append(SNR)
    print('\rlambda = %.2f done'%(lam), end='')

SNR_T = np.sqrt(np.real(inner_p(hS, fr, hS, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)))
snr_INT = snrs/SNR_T

#amps WO
amps = pd.DataFrame(columns=lambdas.round(2))
for n,lam in enumerate(lambdas.round(2)):
    a = []
    for f in fr:
        a.append(AF_PM(f, Ml_WO, zL, y, lam))
    a[0] = 0.
    amps[lam] = a
    print('\r%.1f %% done'%(n/len(lambdas_i)*100), end='')

f_min = 11
f_max = vf_fin(f_min, MTOT_WO/(1+q), MTOT_WO*q/(q+1))

print('\nf max = %.2f'%(f_max))

snrs = []
hS = hp_t_INT * np.conj(amps[1.])
for lam in lambdas_i.round(3):
    hT = hp_t_INT * np.conj(amps[lam])#make_hpL(amps[lam], hp_t_INT, fd=True)
    snr_n = inner_p(hS, fr, hT, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)
    snr_d = inner_p(hT, fr, hT, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)
    SNR = np.real((snr_n/np.sqrt(snr_d)))
    snrs.append(SNR)
    print('\rlambda = %.2f done'%(lam), end='')

SNR_T = np.sqrt(np.real(inner_p(hS, fr, hS, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)))
snr_WO = snrs/SNR_T
#%%
ii = 0
fi = len(lambdas)

# plt.figure()
fig, ax = plt.subplots()#fig.add_subplot(111)

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')

plt.plot(lambdas[ii:fi], np.ones(len(lambdas[ii:fi])), c=dc[3], label='GO')
plt.plot(lambdas[ii:fi], snr_WO[ii:fi], c=dc[2], label='WO')
plt.plot(lambdas_i[ii:], snr_INT[ii:], c=dc[0], label='interference')

plt.xlabel('$\lambda$', fontsize=fontsz)
plt.ylabel('$\\rho/\\rho_{opt}$', fontsize=fontsz)

plt.axhspan(0.998, 1.0005, xmin=0.044, xmax=0.955, color='grey', alpha=0.6, label='threshold')
# plt.axhspan(0.955, 1, xmin=0.046, xmax=0.955, color='limegreen', alpha=0.4, label='th$=0.955$')
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz, labelbottom=True, labeltop=False, top=False)

plt.ylim(0.8, 1.01)
# plt.xlim(0.39, 1.21)
plt.grid()
plt.legend(ncol=4,loc='lower center',bbox_to_anchor=(0.5,1.),framealpha=0.5,)

ii = 51
fi = -11
# '''
inset_axes_1 = inset_axes(ax, width="35%", height="35%", loc=6, axes_kwargs={'alpha':0.5, 'facecolor':'white'})

inset_axes_1.plot(lambdas[ii:fi], np.ones(len(lambdas[ii:fi])), c=dc[3], label='GO')
inset_axes_1.plot(lambdas[ii:fi], snr_WO[ii:fi], c=dc[2], label='WO')
inset_axes_1.plot(lambdas_i[102:-24], snr_INT[102:-24], c=dc[0], label='interference')
plt.xticks([0.92, 1., 1.08])

plt.axhspan(0.998, 1.00005, xmin=0.02, xmax=0.98, color='grey', alpha=0.5, label='threshold') # dc[1]
plt.ylim(0.9969,1.0002)
plt.grid()
plt.tick_params(axis='both',which='both',direction='out',labelsize=9, 
                labelbottom=False, labeltop=True, labelleft=False, labelright=True, 
                top=True, right=True, bottom=False)
# '''
# ax.indicate_inset_zoom(inset_axes_1)#, connector_lines=True)
# mark_inset(ax, inset_axes_1, loc1=4, loc2=2)

plt.show()
# %% save plot
plt.savefig('/home/paolo/Dropbox/PhD/plots/SNR/SNR_regimes.png',
            dpi=300, format='png', bbox_inches="tight") #transparent=True, 
plt.close()












