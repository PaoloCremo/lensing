#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 16:59:35 2021

@author: paolo
"""

#%%load packages and functions
from Dropbox.PhD.Python.program_py.packages.pyCBC_function import *
#%% define case
case = 'd'#'d'
df = 0.2
zS = 0.5
DS = cosmo.luminosity_distance(zS).value #distance of the source in Mpc
zL = 0.01
Ml = 500

MTOT = 100*(1+zS)
q = 1.
m1 = MTOT/(1+q)
m2 = MTOT-m1
hp, hc = get_fd_waveform(approximant='IMRPhenomHM', mass1=m1, mass2=m2,
                         delta_f=df, f_lower=10, inclination=np.pi/3,
                         # s1z = 0.7, s2z = 0.2,
                         distance=DS, coa_phase=0.)

fr = hp.sample_frequencies
ht = hp.data*Fp(0.3, 0.4, 1.5)+hc.data*Fx(0.3, 0.4, 1.5)
# %% defina y e lambda
#old zoom out
# yt = np.arange(0.7, 1.01, 0.05)
# Mt = np.arange(100, 191, 15)

#zoom out
yt = np.arange(0.85, 1.301, 0.05)
Mt = np.arange(350, 620.1, 30)
# yt = np.arange(0.94, 1.081, 0.02)
# Mt = np.arange(450, 521, 10)

#zoom in
#0.75
# yt = np.arange(0.75, 0.8201, 0.01)
# Mt = np.arange(850, 886, 5)
#0.9
# yt = np.arange(0.85, 0.9201, 0.01)
# Mt = np.arange(620, 656, 5)
#1
# yt = np.arange(0.93, 1.0801, 0.02)
# Mt = np.arange(460, 531, 10)
#1.1
# yt = np.arange(1.21, 1.3501, 0.02)
# Mt = np.arange(325, 365, 5)

# lt = np.arange(0.5, 1.21, 0.1)
m1s = len(yt)
print(m1s)
m2s = len(Mt)
print(m2s)
#%% make templates amps
frM = fr[fr<700]
# frM = fr[fr<10**3]
temp_mat = np.ndarray(shape=(m1s,m2s,len(frM)), dtype=np.complex_) # define matrix : y - lambda
now = datetime.now()
lT = 1.
for n1, y in enumerate(yt):
# for n1, Mls in enumerate(Mt):
    # for n2, l in enumerate(lt):
    for n2, Mls in enumerate(Mt):
        for n3, f in enumerate(frM):
            # temp_mat[n1][n2][n3] = AF_PM(f, Ml, zL, y, l)
            # temp_mat[n1][n2][n3] = AF_pm_go(f, Ml, zL, y, l, imag='I')
            # temp_mat[n1][n2][n3] = AF_PM(f, Mls, zL, 1., l)
            temp_mat[n1][n2][n3] = AF_PM(f, Mls, zL, y, lT)
        temp_mat[n1][n2][0] = 0
        print('\r%.f %% done'%((n2+1)/m2s*100), end='')
    print('\r          | y = %.3f done - %.f %% done'%(y, (n1+1)/m1s*100), end='')
    # print('%.1f %% done'%((n1+1)/m1s*100))
    # print('M = %.2f done'%(Mls))

print('\rTemplates done in\n%s'%(datetime.now()-now))

#% make SNR
now = datetime.now()
f_min = 11
f_max = vf_fin(f_min, m1, m2)
print('f max = %.2f'%(f_max))
# Sn = pd.read_csv('/home/paolo/Desktop/AplusDesign.txt', sep='  ', header=None)
dd = 'o3_l1'#'o3_l1'#'aligo'#'aplus'
Sn = pd.read_csv('/home/paolo/Downloads/curves_Jan_2020/'+dd+'.txt', sep=' ', header=None)
Sn.columns = ['fr', 'cs']

snr_mat = np.ndarray(shape=(m1s,m2s))

lam_sig = 1.

amps_p = np.array([AF_PM(f, Ml, zL, 1., lam_sig) for f in fr[fr<700]])
amps_p[0] = 0
hS = make_hpL(amps_p, ht, fd=True)
'''
hS, hc = get_fd_waveform(approximant='IMRPhenomD', mass1=45*(1.1), mass2=45*(1.1),
                             delta_f=0.2, f_lower=10, inclination=0.,
                             distance=D_l(zS, Mpc=True), coa_phase=0.)
hS.resize(41566)
fS = hS.sample_frequencies
hS = hS.data
# '''

pwr = 2
for n1 in range(m1s) : 
    for n2 in range(m2s) :        
        hT = make_hpL(temp_mat[n1][n2], ht, fd=True)
        snr_n = inner_p(hS, fr, hT, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)
        # snr_n = inner_p(hS, fS, hT, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)
        snr_d = inner_p(hT, fr, hT, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)
        '''
        SNR_T = np.sqrt(inner_p(hS, fr, hS, fr, Sn.cs**2, Sn.fr, f_min, f_max))
        snr_mat[n1][n2] = (snr_n/np.sqrt(snr_d))/SNR_T
        '''
        snr_mat[n1][n2] = np.real((snr_n/np.sqrt(snr_d)))
    # print('%.f %% done'%((n1+1)/m1s*100))
    print('\r%.f %% done'%((n1+1)/m1s*100), end='')
    
SNR_T = np.sqrt(np.real(inner_p(hS, fr, hS, fr, Sn.cs**pwr, Sn.fr, f_min, f_max)))
# SNR_T = np.sqrt(inner_p(hS, fS, hS, fS, Sn.cs**pwr, Sn.fr, f_min, f_max))

print('\rSNR done in\n%s'%(datetime.now()-now))
# %% plot
fn = 1
plt.figure() #10x10
# plt.figure(figsize=(8,8))
# plt.figure(figsize=(9,9))
# plt.figure(figsize=(11,11))

plt.matshow(snr_mat/SNR_T, cmap = plt.cm.coolwarm, fignum=fn, vmin=0., vmax=1.)

plt.ylabel('$y$',  fontsize=fontsz)
plt.xlabel('$M_L~[10^2~M_\odot]$',  fontsize=fontsz)

for (i, j), z in np.ndenumerate(snr_mat/SNR_T):
    plt.text(j, i, '{:0.3f}'.format(z), ha='center', va='center', c='k', fontsize=8)
    
cbar = plt.colorbar()#s_m)
cbar.set_label(r'$\rho/\rho_{opt}$', rotation=270,labelpad=15,fontsize=fontsz)
cbar.ax.tick_params(direction='in',labelsize=fontsz)

sg = np.sqrt(11.8)#3
th = 1-(1/2*(sg/SNR_T)**2)

# th = .7
mat = snr_mat/SNR_T
for i in range(m1s):
    for j in range(m2s):
        if mat[i][j] < th:
            mat[i][j] = 0.
        else:
            mat[i][j] = 1.
          
plt.matshow(mat, cmap = plt.cm.Reds, fignum=fn, alpha=0.3) 
plt.xticks(np.arange(m2s), labels=[round(M/100,2) for M in Mt], fontsize=fontsz, rotation=60)
plt.yticks(np.arange(m1s), labels=[round(yy, 3) for yy in yt], fontsize=fontsz)

plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz, labelbottom=True, labeltop=False, top=False)

# plt.title('S/N$_{signal} = %.2f$ - S/N$_{max}= %.2f$'%(SNR_T, snr_mat.max()))
# plt.title('%s - S/N$_{signal} = %.2f$ - %i$\sigma : %.4f$'%(dd,SNR_T, sg, th))
# plt.title(lT)          

plt.show()
print('\nSNR_signal = %.4f\nSNR_max    = %.4f\nthreshold  = %.4f\n'%(SNR_T, snr_mat.max(), 1-(1/2*(sg/SNR_T)**2)))
# %% save
# plt.savefig('/home/paolo/Desktop/waveform/new/Mat_SNR_M'+str(Ml)+'_unl_SM90'+dd+'.png', 
# plt.savefig('/home/paolo/Desktop/waveform/new/Mat_SNR_M'+str(Ml)+'_y'+yst+Temp+'_'+dd+'.png', 
# plt.savefig('/home/paolo/Desktop/waveform/new/Mat_SNR_M'+str(Ml)+'_Ms'+mst+Temp+'_'+dd+'.png', 
# plt.savefig('/home/paolo/Desktop/waveform/new/Mat_SNR_M'+str(Ml)+'_M-y_Ms'+
lTs = str(lT).replace('.', '')
plt.savefig('/home/paolo/Dropbox/PhD/plots/correct/Mat_SNR_M-y_'+lTs+'_Ms500'+
            # mst+Temp+'_L'+str(lam_sig).replace('.','')+'_T075_'+dd+'.png',
            '_q1_'+dd+'_zS05_3sig_2.png',
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()



#%% define parameters
y1 = 1
M1 = 500
l1 = 1

ys = 1
Ms = 500
ls = 0.75

y2 = 4.5
M2 = 250
l2 = 1.1

y3 = 1.63
M3 = 225
l3 = 1

lab1 = '$M_L = %i$\n$y = %.f$\n$\lambda$ = %.f'%(M1, y1, l1)
labs = '$M_L = %i$\n$y = %.f$\n$\lambda$ = %.2f'%(Ms, ys, ls)
lab2 = '$M_L = %i$\n$y = %.2f$\n$\lambda$ = %.f'%(M2, y2, l2)
lab3 = '$M_L = %i$\n$y = %.2f$\n$\lambda$ = %.f'%(M3, y3, l3)
print(lab1)
print('\n'+labs)
print('\n'+lab2)
print('\n'+lab3)
#%% make amps and lensing
amps1 = np.array([AF_PM(f, M1, zL, y1, l1) for f in fr[fr<700]])
ampss = np.array([AF_PM(f, Ms, zL, ys, ls) for f in fr[fr<700]])
amps2 = np.array([AF_PM(f, M2, zL, y2, l2) for f in fr[fr<700]])
amps3 = np.array([AF_PM(f, M3, zL, y3, l3) for f in fr[fr<700]])
amps1[0] = 0
ampss[0] = 0
amps2[0] = 0
amps3[0] = 0
hS1 = make_hpL(amps1, ht, fd=False)
hSs = make_hpL(ampss, ht, fd=False)
hS2 = make_hpL(amps2, ht, fd=False)
hS3 = make_hpL(amps3, ht, fd=False)
dt=10**-4
t = np.arange(0, len(hS1)*dt, dt)
hS1T = make_hpL(amps1, ht, fd=True)
hSsT = make_hpL(ampss, ht, fd=True)
hS2T = make_hpL(amps2, ht, fd=True)
hS3T = make_hpL(amps3, ht, fd=True)
#%% plot TD
norm = 10**-21
plt.figure()
plt.plot(t, hSs/norm, c='k', label=labs)
plt.plot(t, hS1/norm, c='grey', label=lab1)
plt.plot(t, hS2/norm, c='red', label=lab2)
plt.plot(t, hS3/norm, c='b', label=lab3)
# plt.xlim(3.8, 4.03)
plt.grid()
plt.xlabel('t [s]', fontsize=fontsz)
plt.ylabel('strain [$10^{-21}$]', fontsize=fontsz)
plt.legend(ncol=4,loc='lower center',bbox_to_anchor=(0.5,1.),framealpha=0.5, fontsize=fontsz)
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)#, labelbottom=True, labeltop=False)
plt.show()

#%% plot FD

plt.figure()
plt.plot(fr, 2*fr*abs(hSsT), c='k', label=labs)
plt.plot(fr, 2*fr*abs(hS1T), c='grey', label=lab1)
plt.plot(fr, 2*fr*abs(hS2T), c='red', label=lab2)
plt.plot(fr, 2*fr*abs(hS3T), c='b', label=lab3)

plt.xlim(9, 364)
plt.ylim(3*10**-21, 10**-16)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('f [Hz]', fontsize=fontsz)
plt.ylabel('characteristic strain $[\sqrt{4f^2|\\tilde{h}(f)|^2}]$', fontsize=fontsz)
plt.grid()
plt.tick_params(axis='both',which='both',direction='out',labelsize=fontsz)#, labelbottom=True, labeltop=False)
plt.legend()
plt.show()
#%%
# plt.savefig('/home/paolo/Dropbox/PhD/plots/L075_td_1.png',
plt.savefig('/home/paolo/Dropbox/PhD/plots/L075_fd_1.png',
            dpi=300, format='png', transparent=True, bbox_inches="tight")
plt.close()
#%%
