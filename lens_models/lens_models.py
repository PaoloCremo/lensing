#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 19 15:50:03 2020

@author: paolo
"""

#%%load packages and functions
from Dropbox.PhD.Python.program_py.pyCBC_function import *

'''
ATTENZIONE! HERNQUIST ORA HA Mt!!
'''
#%% make lensing
directory = '/home/paolo/Desktop/waveform/'
case = 'h'
lens = 'SIS'
lens_2 = 'NFW'
y = '01'
u1 = '03'
u2 = '-045'
u3 = '-015'
dt = 1
mass = '10-9'
# mass = '5_10-8'     

hp = np.fromfile(directory+'case_'+case+'/hp_td_H74_mt')

file = '/amps_case_h_H74_y'+y+'_M'+mass+'_zL05' # normale
'''
file = '/amps_case_h_H74_y'+y+'_M'+mass+'_zL05_X'+u1+'_'+u2+'_'+u3 # MoG
file_2 = '/amps_case_h_H74_y'+y+'_M'+mass+'_zL05' # normale 2
# '''
if lens == 'Hernquist':
    print('HERN!')
    # file_fd = directory+'case_'+case+'/hp_fd_H74_Mt'     # 32  kB
    file_fd = directory+'case_'+case+'/hp_fd_H74_MT'     # 160 kB
else:
    print('not HENR')
    file_fd = directory+'case_'+case+'/hp_fd_H74_mt'
    
hp_L = make_hpL(directory+'case_'+case+'/'+lens+'/amps/'+file, 
                file_fd)
print('L 1 done')
'''
hp_L_2 = make_hpL(directory+'case_'+case+'/'+lens_2+'/amps/'+file_2, 
                file_fd)
print('L 2 done')
# '''
tU = np.arange(0, len(hp)*dt, dt)

#% roll
pu = np.where(hp==max(hp))[0][0]
pl = np.where(hp_L==max(hp_L))[0][0]
'''
peaks_u = find_peaks(hp)[0]
pl = np.where(hp_L==max(hp_L))[0][0]
peaks_l = find_peaks(hp_L)[0]
hp_l = np.roll(hp_L, pu-pl)
'''
tL =  np.arange(0, len(hp_L)*dt, dt)
# %% plot PLOT plot PLOT
plt.figure()

ee = int(1.94*10**6)
plt.plot(t[:ee], hp[:ee], c='grey', label='unlensed')

# plot_l(lens_2, hp_L_2, 'si', color='k', n_roll = -1302097)
plot_l(lens, tL, hp_L, 'si', pu=pu, pl=pl)#'custom', n_roll = -1042107)

xll, xul = 1.4*10**6, 2.*10**6 
# xll, xul = 0.5*10**6, 7.15*10**6 
plt.xlim(xll, xul) 

ylims = 3.*10**-16
# ylims = 2.*10**-15
# plt.ylim(-ylims, ylims) 

plt.xlabel('t [s]')
plt.ylabel('strain')
plt.grid()
plt.legend(loc='lower left')


'''
print(directory+'case_'+case+'/'+lens+'/plots/M'+mass+'_y'+y+'.png')
# plt.savefig(directory+'case_'+case+'/'+lens+'/plots/M'+mass+'_y'+y+'.png', dpi=300, format='png', bbox_inches="tight")
plt.savefig(directory+'case_'+case+'/'+lens+'/plots/M'+mass+'_y'+y+'_X'+u1+'_'+u2+'_'+u3+'.png', dpi=300, format='png', bbox_inches="tight")

# print(directory+'lens_model/plots/ySIS_01_M'+mass+'/y'+y+'_'+lens+'.png')
# plt.savefig(directory+'lens_model/plots/ySIS_0_M'+mass+'/y'+y+'_'+lens+'.png', dpi=300, format='png', bbox_inches="tight")
'''

plt.show()
# %% save
plt.savefig(directory+'case_'+case+'/'+lens+'/plots/M'+mass+'_y'+y+'.png', dpi=300, format='png', bbox_inches="tight")
# plt.savefig(directory+'case_'+case+'/'+lens+'/plots/M'+mass+'_y'+y+'_X'+u1+'_'+u2+'_'+u3+'.png', dpi=300, format='png', bbox_inches="tight")
# plt.savefig(directory+'lens_model/plots/ySIS_1_M'+mass+'/y'+y+'_'+lens+'_p.png', dpi=300, format='png', bbox_inches="tight")

#%% find maxima & zeros
# pL = find_peaks(abs(np.roll(hp_L, pu-pl)))[0]
hp_l = np.roll(hp_L, pu-pl)
pL = find_peaks(abs(hp_l))[0]
pU_p = find_peaks(abs(hp))[0]
pU = []
for n,i in enumerate(pU_p[:-1]):
    # print(t[pU_p[n+1]]-t[i])
    if (tU[pU_p[n+1]]-tU[i])>5000:
        pU.append(i)

zU = find_zeros(hp)
zL = find_zeros(hp_l)
tp = tU[pu]
zUm = np.where(zU>tp)[0][0]
zLm = np.where(zL>tp)[0][0]
#%% find max(peaks)
mu = np.where(abs(hp.data) == max(abs(hp)))[0][0]
mu = np.where(pU == mu)[0][0]
ml = np.where(abs(hp_l) == max(abs(hp_l)))[0][0]
ml = np.where(pL == ml)[0][0]
#%% phase shift plot
# ''' # zeros
plt.figure()
y_pl = (tU[zU[-40+zUm:zUm+5]]-tL[zL[-40+zLm:zLm+5]])/tU[zU[-40+zUm:zUm+5]]*100

plt.scatter(tU[zU[-40+zUm:zUm+5]], y_pl)

plt.vlines(tU[zU[zUm]], min(y_pl), max(y_pl), label='merger')
plt.xlim(tU[zU[0]], tU[zU[zUm+8]])
''' # max
plt.figure()
y_pl = (t[pU[-40+mu:mu+5]]-tL[pL[-40+ml:ml+5]])/t[pU[-40+mu:mu+5]]*100

plt.scatter(t[pU[-40+mu:mu+5]], y_pl)
plt.vlines(t[pU[mu]], min(y_pl), max(y_pl), label='merger')
plt.xlim(t[pU[0]], t[pU[mu+8]])
# '''

plt.xlabel('peak time unlensed [s]')
plt.ylabel('% difference (zeros)')#' from no phase shift')
# plt.ylabel('ratio')#' from no phase shift')
# plt.xscale('log')
plt.grid()
plt.legend()
plt.show()
#%%
# plt.savefig(directory+'lens_model/plots/features/ps0_'+lens+'_M'+mass+'_y'+y+'_z.png', dpi=300, format='png', bbox_inches="tight")
plt.savefig(directory+'lens_model/plots/deg/ps0_feat_'+ss+'_M'+mass+'_y'+y+'_z02.png', dpi=300, format='png', bbox_inches="tight")

#%% peak ratios plot
y_to_pl =  abs(hp[pU[-40+mu:mu+5]])/abs(hp_l[pL[-40+ml:ml+5]])
# y_to_pl =  abs(hp_l[pL[-40+ml:ml+5]])/abs(hp[pU[-40+mu:mu+5]])
plt.figure()
'''
# plt.plot(range(45), abs(hp[pU[-40+mu:mu+5]]*6), c='grey', label='unlensed')
# plt.plot(range(45), abs(hp_l[pL[-40+ml:ml+5]]), c='orangered', label='lensed')
# plt.xlabel('peak count')
# plt.ylabel('peak strain')
# '''
# plt.plot(range(43), abs(hp[pU[-40+mu:mu+3]])/abs(hp_l[pL[-40+ml:ml+3]]))
plt.plot(tU[pU[-40+mu:mu+5]], y_to_pl)

plt.vlines(tU[pU[mu]], min(y_to_pl), max(y_to_pl), label='merger')
plt.xlabel('time unlensed [s]')
plt.ylabel('peak ratio')
# '''
plt.grid()
plt.legend()#loc='lower left')
plt.show()
#%%
# plt.savefig(directory+'lens_model/plots/features/pr_'+lens+'_M'+mass+'_y'+y+'_2.png', dpi=300, format='png', bbox_inches="tight")
# plt.savefig(directory+'lens_model/plots/features/pstr_'+lens+'_M'+mass+'_y'+y+'.png', dpi=300, format='png', bbox_inches="tight")
plt.savefig(directory+'lens_model/plots/deg/pr_feat_'+ss+'_M'+mass+'_y'+y+'_z02.png', dpi=300, format='png', bbox_inches="tight")

# %%
#  M=10^9
rSIS = 9.62299*10**18 # theta_E * Dl
rNIS = 9.62299*10**18
rNFW = 9.065935*10**20
rNF2 = 3.4313*10**20
rHER = 3.5948*10**20
rGN0 = 3.3531*10**20
rGN2 = 1.517823*10**20
rBUR = 3.43436*10**20
#%%
#  M=5*10^8
rSIS = 6.804480277940293*10**18 #d
rNIS = 6.804480277940293*10**18 #d
rNFW = 8.39304*10**20
rNF2 = 2.88202*10**20
rHER = 3.32636*10**20
rGN0 = 2.82031*10**20
rGN2 = 1.2713*10**20
rBUR = 2.88819421610*10**20
#%%
rs = [rSIS, rNIS, rNFW, rNF2, rHER, rGN0, rGN2, rBUR]
rs_n = ['ySIS', 'yNIS', 'yNFW', 'yNF2', 'yHER', 'yGN0', 'yGN2', 'yBUR']
#%
y_SIS = 1.
ys = [y_SIS*rs[0]/r for r in rs]
[print('%s = %.5f'%(n, np.round(yy, 5))) for n,yy in zip(rs_n,ys)] 
#%% frequency domain
par = pd.read_csv('/home/paolo/Desktop/waveform/lens_model/codes/parameters.txt', index_col=False)
# par = pd.read_csv('/home/paolo/Desktop/waveform/lens_model/codes/parameters_MoG.txt', index_col=False)
directory = '/home/paolo/Desktop/waveform/'
case = 'h'
# lens = 'Burkert'
sce = 7
y = par.SIS[sce]
dt = 1
if sce % 2 == 0 :
    mass = '10-9'
else : 
    mass = '5_10-8'

print('M = '+ mass + '\ny = ' + y[:1] + '.' + y[1:])
#%%
hp_t = np.fromfile(directory+'case_'+case+'/hp_fd_H74_mt')
l = len(hp_t)-2
hp_t = unite(hp_t)
fr = np.arange(l/2+1)/(1*l)

#%%plot
ii = 29
ff = 635
plt.figure()#figsize=(12.8, 9.6))

plt.plot(fr[ii:ff], 2*fr[ii:ff]*abs(hp_t[ii:ff]), c='grey', linestyle='solid',label='unlensed')
#%
for ss in par.columns[4:]:
    file = '/amps_case_h_H74_y'+par[ss][sce]+'_M'+mass+'_zL05'
    if ss == 'Hernquist':
        #  scenario 6/7
        # file_hp = directory+'case_'+case+'/hp_fd_H74_Mt' #  32 kB
        # lfs = int(10**7)
        # ii = 95
        # ff = 2005
        file_hp = directory+'case_'+case+'/hp_fd_H74_MT' #160 kB
        lfs = int(5*10**7)
        ii = 450
        ff = 9975
        fs = np.arange(lfs/2+1)/(1*lfs)
    else : 
        file_hp = directory+'case_'+case+'/hp_fd_H74_mt'
        fs = fr
        ii = 29
        ff = 635
    hp_L_t = make_hpL(directory+'case_'+case+'/'+ss+'/amps/'+file, 
                file_hp, fd = True)
    lab = ss 
    
    plt.plot(fs[ii:ff], 2*fr[ii:ff]*abs(hp_L_t[ii:ff]), c = par[ss][0], 
             linestyle = par[ss][1], label=lab)
    print(ss + ' DONE!')

plt.legend(fontsize=8)#, loc='upper left')
#%
plt.xlim(10**-5, 1.1*10**-4)

# plt.ylim(1.5*10**-16, 10**-15) # y = 1
# plt.ylim(1.5*10**-16, 10**-14)   # y = 0.1
# '''
# plt.ylim(1.5*10**-16, 10**-13) # y = 1
# plt.ylim(1.5*10**-16, 2*10**-13) # y = 0.1
plt.ylim(1.5*10**-16, 3*10**-13) # y = 0
# '''
plt.xscale('log')
plt.yscale('log')
plt.xlabel('f [Hz]')
plt.ylabel('characteristic strain $[\sqrt{4f^2|\\tilde{h}(f)|^2}]$')
# handles, labels = plt.gca().get_legend_handles_labels()
# order = [0, 2, 3, 4, 5, 6, 7, 1, 8]
# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=8)
# major_t = [x*10**-5 for x in range(10)]
plt.grid(True, which='both', ls='-')
# plt.set_xdata=major_t

'''
plt.savefig(directory+'lens_model/plots/ySIS_'+par.SIS[sce]+'_M'+mass+'/fd_1.png', dpi=300, format='png', bbox_inches="tight")
plt.savefig(directory+'lens_model/plots/ySIS_'+par.SIS[sce]+'_M'+mass+'/fd_2.png', dpi=300, format='png', bbox_inches="tight")
'''
# plt.savefig(directory+'case_h/Hernquist/plots/y_'+ys[0]+'_M'+mass+'_fd.png', dpi=300, format='png', bbox_inches="tight")

plt.show()

#%% frequency domain MoG
par = pd.read_csv('/home/paolo/Desktop/waveform/lens_model/codes/parameters_MoG.txt', index_col=False)
directory = '/home/paolo/Desktop/waveform/'
case = 'h'
# lens = 'Burkert'
sce = 7
y = par.SIS[sce]
dt = 1
if sce % 2 == 0 :
    mass = '10-9'
else : 
    mass = '5_10-8'

print('M = '+ mass + '\ny = ' + y[:1] + '.' + y[1:])
#%%
hp_t = np.fromfile(directory+'case_'+case+'/hp_fd_H74_mt')
l = len(hp_t)-2
hp_t = unite(hp_t)
fr = np.arange(l/2+1)/(1*l)
#%%plot MoG
ii = 29
ff = 635
plt.figure()#figsize=(12.8, 9.6))

plt.plot(fr[ii:ff], 2*fr[ii:ff]*abs(hp_t[ii:ff]), c='grey', linestyle='solid',label='unlensed')
#%
for ss in par.columns[1:]:
    
    if 'p' in ss :
        file = '/amps_case_h_H74_y'+par[ss][sce]+'_M'+mass+'_zL05_X-03_045_015'
        mod = ss[:-2]
    elif 'm' in ss :
        file = '/amps_case_h_H74_y'+par[ss][sce]+'_M'+mass+'_zL05_X03_-045_-015'
        mod = ss[:-2]
    else : 
        file = '/amps_case_h_H74_y'+par[ss][sce]+'_M'+mass+'_zL05'
        mod = ss        
   
    file_hp = directory+'case_'+case+'/hp_fd_H74_mt'
    hp_L_t = make_hpL(directory+'case_'+case+'/'+mod+'/amps/'+file, 
                file_hp, fd = True)
    lab = ss 
    fs = fr
    plt.plot(fs[ii:ff], 2*fr[ii:ff]*abs(hp_L_t[ii:ff]), c = par[ss][0], 
             linestyle = par[ss][1], label=lab)
    print(ss + ' DONE!')

plt.legend(fontsize=6)#, loc='lower left')
#%
plt.xlim(10**-5, 1.1*10**-4)

plt.ylim(1.5*10**-16, 2.*10**-14) # y = 0
# '''
plt.xscale('log')
plt.yscale('log')
plt.xlabel('f [Hz]')
plt.ylabel('characteristic strain $[\sqrt{4f^2|\\tilde{h}(f)|^2}]$')
# handles, labels = plt.gca().get_legend_handles_labels()
# order = [0, 2, 3, 4, 5, 6, 7, 1, 8]
# plt.legend([handles[idx] for idx in order],[labels[idx] for idx in order], fontsize=8)
# major_t = [x*10**-5 for x in range(10)]
plt.grid(True, which='both', ls='-')
# plt.set_xdata=major_t

'''
plt.savefig(directory+'lens_model/plots/ySIS_'+par.SIS[sce]+'_M'+mass+'/fd_MoG.png', dpi=300, format='png', bbox_inches="tight")
'''
# plt.savefig(directory+'case_h/Hernquist/plots/y_'+ys[0]+'_M'+mass+'_fd.png', dpi=300, format='png', bbox_inches="tight")

plt.show()


#%% different ys
directory = '/home/paolo/Desktop/waveform/'
case = 'h'
lens = 'MoG_2'
# sce = 7
# ys = ['0', '01', '1', '5']
ys = ['08', '09', '1', '11', '12']
u1 = '-03'
u2 = '045'
u3 = '015'
dt = 1
mass = '10-9'
# mass = '5_10-8'

print('M = '+ mass + '\nlens : ' + lens)
#%%
hp_t = np.fromfile(directory+'case_'+case+'/hp_fd_H74_mt')
l = len(hp_t)-2
hp_t = unite(hp_t)
fr = np.arange(l/2+1)/(1*l)
#%%% plot
ii = 29
ff = 635
plt.figure()#figsize=(12.8, 9.6))
plt.title(lens)
plt.plot(fr[ii:ff], 2*fr[ii:ff]*abs(hp_t[ii:ff]), c='grey', linestyle='solid',label='unlensed')

for y in ys:
    # file = '/amps_case_h_H74_y'+y+'_M'+mass+'_zL05'
    file = '/amps_case_h_H74_y'+y+'_M'+mass+'_zL05_X'+u1+'_'+u2+'_'+u3 # MoG
    print('\n' + file)
    # ''' # normal
    file_hp = directory+'case_'+case+'/hp_fd_H74_mt'
    fs = fr
    ii = 29
    ff = 635
    ''' # Hernquist
    if y == '01':
        file_hp = directory+'case_'+case+'/hp_fd_H74_MT'
        lfs = int(5*10**7)
        ii = 450
        ff = 9975
    else:
        file_hp = directory+'case_'+case+'/hp_fd_H74_Mt'
        lfs = int(10**7)
        ii = 95
        ff = 2005
        
    fs = np.arange(lfs/2+1)/(1*lfs)
    '''
    hp_L_t = make_hpL(directory+'case_h/'+lens+'/amps/'+file, 
                file_hp, fd = True)
    if len(y)==1:
        lab = 'y = ' + y
    else :
        lab = 'y = ' + y[:1] + '.' + y[1:]
    
    plt.plot(fs[ii:ff], 2*fr[ii:ff]*abs(hp_L_t[ii:ff]), label=lab)
    print(y + ' DONE!')

plt.legend(fontsize=8)#, loc='upper left')

plt.xlim(10**-5, 1.1*10**-4)
plt.ylim(1.5*10**-16, 4*10**-16) # 
# plt.ylim(1.5*10**-16, 6*10**-15) # main
'''
plt.ylim(1.5*10**-16, 5*10**-15) # SIS  
plt.ylim(2.*10**-16, 3.7*10**-16) # NIS & GNFW_0 & Burkert 
plt.ylim(1.5*10**-16, 6*10**-15) # NFW & GNFW_2 & main
plt.ylim(1.7*10**-16, 5*10**-16) # NFW_2 
plt.ylim(1.5*10**-16, 3*10**-13) # Hernquist
plt.ylim(1.*10**-17, 2*10**-14) # MoG

plt.legend(fontsize=8)
'''
plt.xscale('log')
plt.yscale('log')
plt.xlabel('f [Hz]')
plt.ylabel('characteristic strain [$\sqrt{4f^2|\\tilde{h}(f)|^2}$]')

# https://arxiv.org/pdf/1408.0740.pdf

plt.grid(True, which='both', ls='-')
'''
MAIN
plt.savefig(directory+'case_h/fd_plots/'+lens+'_M'+mass+'_fd.png', dpi=300, format='png', bbox_inches="tight")
ZOOM
plt.savefig(directory+'case_h/fd_plots/'+lens+'_M'+mass+'_fd_z.png', dpi=300, format='png', bbox_inches="tight")
MOG
plt.savefig(directory+'case_h/fd_plots/'+lens+'_M'+mass+'_fd_z_X'+u1+'_'+u2+'_'+u3+'.png', dpi=300, format='png', bbox_inches="tight")

plt.savefig(directory+'case_'+case+'/'+lens+'/plots/fd/M'+mass+'_fd_U'+u1+'_'+u2+'_'+u3+'.png', dpi=300, format='png', bbox_inches="tight")

'''
plt.show()
#%%
