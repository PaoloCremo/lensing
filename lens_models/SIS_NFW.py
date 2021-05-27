#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  1 18:03:56 2020

@author: paolo
"""

#%%load packages and functions
from Dropbox.PhD.Python.program_py.packages.pyCBC_function import *

def subplot(axs, i, j, file_name, color, lab, lw=1, roll=True):
    
    dt = 1
    
    if not 'lensed' in file_name:
        amps = np.fromfile(file_name, np.complex64)
        if len(amps) == 2157:
            hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74'
            hp_t_fn = '/home/paolo/Desktop/waveform/case_h/hp_fd_H74'
            l = 2157032
        elif len(amps) == 1262:
            hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74_mt'
            hp_t_fn = '/home/paolo/Desktop/waveform/case_h/hp_fd_H74_mt'
            l = 3157032
                
        hp_t = np.fromfile(hp_t_fn)
        hp_t = unite(hp_t)
        amps = np.concatenate((amps, np.zeros(len(hp_t)-len(amps))))
        hp_L = make_transform_with_C_c2r(hp_t * np.conj(amps))
        hp_L /= len(hp_L)
    else:
        hp_L = np.fromfile(file_name)
        hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74_mt'
    # np.roll(hp_L, shift)
    hp = np.fromfile(hp_fn)
    t = np.arange(0, len(hp)*dt, dt)
    pu = np.where(hp==max(hp))[0][0]
    pl = np.where(hp_L==max(hp_L))[0][0]
    if roll:
        hp_L = np.roll(hp_L, int(pu-pl))
    if file_name[36:39] == 'NFW':
        if lab == 'y = 7':
            hp_L = np.roll(hp_L, int(pu-2.3842*10**6))
        if lab == 'y = 0.1':
            hp_L = np.roll(hp_L, int(pu -1.86607*10**6))
    axs[i,j].plot(t, hp, c='grey', label='unlensed') 
    axs[i,j].plot(t, hp_L, linewidth=lw, c=color, label=lab)
    axs[i,j].grid()
    axs[i,j].legend(loc='upper left')
    
def subplot_fd(axs, i, j, file_name, color, lab, u_l=10**-14, lw=1):
    
    amps = np.fromfile(file_name, np.complex64)
    if len(amps) == 2157:
        hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74'
        hp_t_fn = '/home/paolo/Desktop/waveform/case_h/hp_fd_H74'
        l = 2157032
    elif len(amps) == 1262:
        hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74_mt'
        hp_t_fn = '/home/paolo/Desktop/waveform/case_h/hp_fd_H74_mt'
        l = 3157032
        
    dt = 1
    hp_t = np.fromfile(hp_t_fn)
    hp_t = unite(hp_t)
    amps = np.concatenate((amps, np.zeros(len(hp_t)-len(amps))))
    hp_L_t = hp_t * np.conj(amps)
    fr = np.arange(l/2+1)/(dt*l)
    '''
    hp_L = make_transform_with_C_c2r(hp_t * np.conj(amps))
    hp_L /= len(hp_L)
    # np.roll(hp_L, shift)
    hp = np.fromfile(hp_fn)
    t = np.arange(0, len(hp)*dt, dt)
    pu = np.where(hp==max(hp))[0][0]
    pl = np.where(hp_L==max(hp_L))[0][0]
    if pl < l/2:
        hp_L = np.roll(hp_L, int(3*l/4-pl))
        pl = np.where(hp_L==max(hp_L))[0][0]
    # if lab == 'y = 0.1':
        # pl += (1.86607-1.85113)*10**6
    '''
    axs[i,j].plot(fr, 2*fr*abs(hp_t), c='grey', label='unlensed') 
    axs[i,j].plot(fr, 2*fr*abs(hp_L_t), linewidth=lw, c=color, label=lab)
    axs[i,j].set_xscale('log')
    axs[i,j].set_yscale('log')
    axs[i,j].set_ylim(10**-20, u_l)
    axs[i,j].grid()
    axs[i,j].legend()

def make_hpL(fn, fn_ut=0, fd=False):
    amps = np.fromfile(fn, np.complex64)
    if len(amps) == 2157:
        hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74'
        hp_t_fn = '/home/paolo/Desktop/waveform/case_h/hp_fd_H74'
        l = 2157032
    elif len(amps) == 1262:
        hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74_mt'
        hp_t_fn = '/home/paolo/Desktop/waveform/case_h/hp_fd_H74_mt'
        l = 3157032
    hp_t = np.fromfile(fn_ut)
    hp_t = unite(hp_t)
    amps = np.concatenate((amps, np.zeros(len(hp_t)-len(amps))))
    hp_L_t = hp_t * np.conj(amps)
    hp_L = make_transform_with_C_c2r(hp_L_t)
    hp_L /= len(hp_L)
    
    if not fd:
        return hp_L
    else:
        return hp_L_t

def plot_vs(fn_SIS, fn_NFW, ylab1, ylab2, save_file = False):
    
    hp = np.fromfile('/home/paolo/Desktop/waveform/case_h/hp_td_H74')
    hp_L_SIS = make_hpL(fn_SIS)
    hp_L_NFW = make_hpL(fn_NFW)
    dt = 1
    t = np.arange(0, len(hp_L_SIS)*dt, dt)
    # tNFW = np.arange(0, len(hp_L_NFW)*dt, dt)
    pu = np.where(hp==max(hp))[0][0]
    pSIS = np.where(hp_L_SIS==max(hp_L_SIS))[0][0]
    pNFW = np.where(hp_L_NFW==max(hp_L_NFW))[0][0]
    # pNFW += (1.86607-1.85113)*10**6 # case 0.1
    hp_L_SIS = np.roll(hp_L_SIS, int(pu-pSIS))
    hp_L_NFW = np.roll(hp_L_NFW, int(pu-pNFW))
    hp_L_NFW = hp_L_NFW[:len(hp_L_SIS)]
    
    
    plt.figure(figsize=(16,9))
    
    plt.plot(t, hp, c='grey', label='grey')
    plt.plot(t, hp_L_SIS, c='forestgreen', label='SIS - y = %s'%(ylab1))
    plt.plot(t, hp_L_NFW, c='orangered', label='NFW - y = %s'%(ylab2))
    
    plt.legend(loc='upper left')
    plt.grid()

    # plt.xlim(1.72*10**6, 1.9*10**6)
    plt.xlim(1.72*10**6, 1.94*10**6)
    
    if save_file:
        plt.savefig(save_file, dpi=300, format='png', bbox_inches="tight")
    plt.show()

def plot_vs_fd(fn_SIS, fn_NFW, ylab1, ylab2, save_file = False):
    
    hp_L_SIS = make_hpL(fn_SIS, True)
    hp_L_NFW = make_hpL(fn_NFW, True)
    hp_t = np.fromfile('/home/paolo/Desktop/waveform/case_h/hp_fd_H74')
    hp_t = unite(hp_t)
    lSIS = 2*len(hp_L_SIS)-2
    lNFW = 2*len(hp_L_NFW)-2
    
    
    dt = 1
    frSIS = np.arange(lSIS/2+1)/(lSIS*dt)
    frNFW = np.arange(lNFW/2+1)/(lNFW*dt)

    plt.figure(figsize=(16,9))
    
    if len(frSIS) == len(hp_t):
        plt.plot(frSIS, 2*frSIS*abs(hp_t), c='grey', label='unlensed')
    else:
        plt.plot(frNFW, 2*frNFW*abs(hp_t), c='grey', label='unlensed')
    plt.plot(frSIS, 2*frSIS*abs(hp_L_SIS), c='forestgreen', label='SIS - y = %s'%(ylab1))
    plt.plot(frNFW, 2*frNFW*abs(hp_L_NFW), c='orangered', label='NFW - y = %s'%(ylab2))
    plt.xscale('log')
    plt.yscale('log')
    plt.grid()
    plt.legend(loc='lower left')
    plt.xlim(8*10**-6, 2*10**-4)
    plt.ylim(10**-20, 10**-14)
    
    
    if save_file:
        plt.savefig(save_file, dpi=300, format='png', bbox_inches="tight")
    plt.show()


#%%
cc = 'orangered'
fig, axs = plt.subplots(4, 2, figsize=(16,20), sharex=True, gridspec_kw={'hspace': 0.15, 'wspace': 0.1})

ys = ['0', '0001', '001', '01', '1', '3', '5', '7']
legends =  ['0', '0.001', '0.01', '0.1', '1', '3', '5', '7']
# for y, ii, jj,ll in zip(ys,[0,0,1,1,2,2],[0,1,0,1,0,1],legends):
for y, ii, jj,ll in zip(ys,[0,1,2,3,0,1,2,3],[0,0,0,0,1,1,1,1],legends):
    # f_n = '/home/paolo/Desktop/waveform/case_h/Hernquist/amps/amps_case_h_H74_y'+y+'_M10-9_zL05'
    f_n = '/home/paolo/Desktop/waveform/case_h/Hernquist/lensed/hpL_case_h_H74_y'+y+'_M10-9_zL05'
    # print('\n', y)
    # print('lensed' in f_n)
    subplot(axs, ii, jj, f_n, cc, 'y = %s'%(ll), roll=y not in ['0001', '001', '01','7'])
    # if not (ii == 3 and jj == 1):
        # axs[ii,jj].set_xlim(1.72*10**6, 1.9*10**6)
    
    # if jj == 0 :
        # axs[ii,jj].set_ylim(-1.7*10**-15, 1.7*10**-15)
    # elif jj == 1:
        # axs[ii,jj].set_ylim(-1.16*10**-15, 1.16*10**-15)

# axs[3,1].set_xlim(1.25*10**6, 2.75*10**6)

# plt.xlim(1.72*10**6, 1.94*10**6) #SIS
# plt.xlim(1.72*10**6, 1.9*10**6) #NFW

# fig.savefig('/home/paolo/Desktop/waveform/lens_model_plots/Hernquist_td.png', dpi=300, format='png', bbox_inches="tight")
#%%
cc = 'orangered'
fig, axs = plt.subplots(4, 2, figsize=(16,20), sharex=True, gridspec_kw={'hspace': 0, 'wspace': 0.1})

ys = ['0', '0001', '001', '01', '1', '3', '5', '7']
legends =  ['0', '0.001', '0.01', '0.1', '1', '3', '5', '7']

for y, ii, jj,ll in zip(ys,[0,1,2,3,0,1,2,3],[0,0,0,0,1,1,1,1],legends):
    f_n = '/home/paolo/Desktop/waveform/case_h/Hernquist/amps/amps_case_h_H74_y'+y+'_M10-9_zL05'
    subplot_fd(axs, ii, jj, f_n, cc, 'y = %s'%(ll), u_l=10**-13)

plt.xlim(8*10**-6, 2*10**-4)

# fig.savefig('/home/paolo/Desktop/waveform/lens_model_plots/Hernquist_fd.png', dpi=300, format='png', bbox_inches="tight")

#%%
ys = ['0', '0001', '001', '01', '1', '10']
legends =  ['0', '0.001', '0.01', '0.1', '1', '10']
i = 5
j = 3
fn_SIS = '/home/paolo/Desktop/waveform/case_h/SIS/amps/amps_case_h_H74_y'+ys[i]+'_M10-9_zL05'
fn_NFW = '/home/paolo/Desktop/waveform/case_h/NFW/amps/amps_case_h_H74_y'+ys[j]+'_M10-9_zL05'
plot_vs_fd(fn_SIS, fn_NFW, legends[i], legends[j])#, '/home/paolo/Desktop/waveform/lens_model_plots/SIS_y'+ys[i]+'vsNFW_y'+ys[j]+'_fd.png')

#%%
zS = 1.
zL = 0.5

Dl = D_l(zL, D_a=True)
Ds = D_l(zS, D_a=True)
Dls = D_l(zS, zL, D_a=True)
#%% 
M = 10**9 # = M_200
# M_200 = 8*M

# NFW
c_200 = 10

# rho_c = 1.8788*10**-26 * H0/100 # 3*H^2/(8pi G)
rho_c = 3*(1/Hz_1(zL, H0)/pctomt/10**3)**2/(8*np.pi*G)
rho_s = 200/3*rho_c*c_200**3/(np.log(1+c_200)-c_200/(1+c_200)) #NFW
rho_s_hern = 2*200/3*rho_c*c_200*(1+c_200) #Hernquist 
rho_S_GNFW = 200/3*rho_c*(c_200**3)/np.log(1+c_200)

def OM(zz):
    om = (Om*(1+zz)**3) + Og*(1+zz)**4 + Ok*(1+zz)**2 + Ol*(1+zz)**(3*(1+w_0+w_a))
    
    return om
                                                                    

def h(x):
    h = -((2*cmath.sqrt(-1 + 2/(1+x)) * cmath.atanh(cmath.sqrt(-1+2/(1+x))))/(-1+x))+np.log(x/2)
    
    return h

h1 = np.real(np.mean([h(0.9999), h(1.0001)]))
rs = (M*smtokg/(4*np.pi*rho_s*h1))**(1/3)

def M_NFW(x, rhos=rho_s, rs=rs):
    M = np.real(4*np.pi*rhos*rs**3*h(x))
    
    return M

# SIS

def sig(M):
    s = ((M*smtokg*c**2*G)/(4*np.pi**2)*(Ds/(Dl*Dls)))**(1/4)
    
    return s

def tE(ss):
    te = 4*np.pi*(ss/c)**2*Dls/Ds
    
    return te
# '''
def tE2(Mass):
    te = np.sqrt(4*G*Mass*smtokg/c**2*Dls/(Dl*Ds))
    
    return te
# '''   
    
def M_SIS_te(ssig):
    m = 4*np.pi**2/G*ssig**4/c**2*Dl*Dls/Ds
    return m
    
def M_SIS(x, m=10**9) :
    M = x*m*smtokg#Dl*tE*2*sigma**2/G
    
    return M


#%%
# rs = 4.14195*10**20
# rs = 2.3720433*10**20
rs = 9.065935*10**20

r = np.concatenate((np.linspace(10**-3, 0.99, 10**4),
                     np.linspace(1.01, 10**2, 10**4)))
y_NFW = [M_NFW(xx, rs=rs)/smtokg for xx in r]
y_SIS = [M_SIS(xx)/smtokg for xx in r]
'''
r = np.concatenate((np.linspace((3*10**17)/pctomt, (0.99*rs)/pctomt, 10**4),
                    np.linspace((1.01*rs)/pctomt, (3*10**22)/pctomt, 10**4)))
y_NFW = [M_NFW(rr*pctomt/rs)/smtokg for rr in r]
y_SIS = [M_SIS(rr*pctomt/(Dl*tE(sig(M))))/smtokg for rr in r]
'''
tE_NFW = tE2(M)
r_te = tE_NFW*Dl/rs
#%%
plt.figure(figsize=(8,4.5))

plt.title('$M = %.e$ - $c_{200} = %.f$ - $\\rho_s = %.2e$ - $r_s = %.2e$\n$M_{NFW}(y=1) = %.2e$ - $\\theta_E*D_l/r_s = %.3e$'%(M, c_200, rho_s, rs, M_NFW(1.0000001, rs=rs)/smtokg, r_te))

lw=3
plt.plot(r, y_NFW, c='b', linewidth=lw, label='NFW\n$r_s$')
plt.plot(r, y_SIS, c='red', linewidth=lw, label='SIS\n$\\theta_E$')

plt.vlines(r_te, 0, max(y_NFW+y_SIS), color='b')#, label='$r_s$')
# plt.vlines(tE(sig(M))*Dl/pctomt, 0, max(y_NFW+y_SIS), color='red')#, label='$\\theta_E$')

plt.hlines(10**9, r[0], r[-1], label='$10^9~M_\odot$')

plt.xscale('log')
plt.yscale('log')

# plt.xlabel('r [pc]')
plt.xlabel('$\\xi/ \\xi_0$')
plt.ylabel('M [M$_\odot$]')

plt.legend()
plt.grid()

Mstring = '{:.1e}'.format(M).replace('.', '_')
file = '/home/paolo/Desktop/waveform/case_h/NFW/plots/SISvsNFW_M'+Mstring+'_c_ce2'+str(int(c_200))+'.png'
print('M = %.e\nc = %.f'%(M, c_200))
print(file)
# '''
plt.savefig(file[:-4]+'_pE2.png', dpi=300, format='png', bbox_inches="tight")
# '''
plt.show()
#%%
tE_NFW = tE2(M)
r_te = tE_NFW*Dl/rs

#%%
# rs = 2.3720433*10**20
# rs = 9.065935*10**20
rs = 9.065935*10**20
r_te = tE_NFW*Dl/rs
print('\n%.5e\n%.6e'%(r_te, M_NFW(r_te, rs=rs)/smtokg))


#%%
'''
dt = 1
hp = np.fromfile('/home/paolo/Desktop/waveform/case_h/hp_td_H74')
t = np.arange(0, len(hp)*dt, dt)
hp_mt = np.fromfile('/home/paolo/Desktop/waveform/case_h/hp_td_H74_mt')
t_mt = np.arange(0, len(hp_mt)*dt, dt)
''




