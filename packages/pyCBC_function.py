#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 15:27:50 2020

@author: paolo
"""
# https://lscsoft.docs.ligo.org/lalsuite/lalsimulation/group___l_a_l_sim_i_m_r_phenom__c.html

import os
import math
import cmath
import pylab
# import quadpy
import numpy as np
import pandas as pd
from mpmath import hyp1f1
# from scipy.special import hyp
from scipy.special import gamma
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.signal import find_peaks
from scipy.interpolate import interp1d
from datetime import datetime, timedelta
from astropy.cosmology import Planck18_arXiv_v2 as cosmo

from pycbc import frame
from pycbc.types import real_same_precision_as
from pycbc.filter import resample_to_delta_t, highpass, matched_filter, sigma
from pycbc.types import TimeSeries, FrequencySeries, Array
from pycbc.waveform import get_td_waveform, get_fd_waveform, td_approximants
from pycbc.psd import interpolate, inverse_spectrum_truncation
from pycbc.waveform.utils import time_from_frequencyseries, phase_from_frequencyseries, phase_from_polarizations

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

# define functions

def Hz_1(zz, H0, Om=Om, Og=Og, Ok=Ok, Ol=Ol, w_0=w_0, w_a=w_a):
    '''
    calculate Hubble parameter ^ -1
    '''
    hz =  1/(math.sqrt((H0**2)*((Om*(1+zz)**3) + Og*(1+zz)**4 + Ok*(1+zz)**2 + Ol*(1+zz)**(3*(1+w_0+w_a))*math.exp(-3*w_a*zz/(1+zz)))))

    return hz

def D_l(zz, z_in=0, Mpc=False, c=c/10**3, H0=H0, Om=Om, Og=Og, Ok=Ok, Ol=Ol, w_0=w_0, w_a=w_a, Hz_1=Hz_1, D_a=False):
    '''
    calculate luminosity distance
    '''
    arg = (zz, H0, Om, Og, Ok, Ol, w_0, w_a)
    dist = ((1+zz)*c)*integrate.quad(Hz_1, z_in, zz, args=(arg[1:]))[0]
    
    if D_a:
        dist = (c/(1+zz))*integrate.quad(Hz_1, z_in, zz, args=(arg[1:]))[0]
    
    if Mpc:
        return dist
    elif not Mpc:
        return dist*3.08567758128*10**22 #in meters

def make_wf(m1, m2, dt, f_i, inc, dist, appr='SEOBNRv4T'):
    hp, hc = get_td_waveform(approximant=appr, mass1=m1, mass2=m2,
                             delta_t=dt, f_lower=f_i, inclination=inc,
                             distance=dist)
    
    return hp, hc

def make_wf_tilde(m1, m2, df, f_i, inc, dist, appr='SEOBNRv4T'):
    hp_tilde, hc_tilde = get_fd_waveform(approximant=appr, mass1=m1, mass2=m2,
                             delta_f=df, f_lower=f_i, inclination=inc,
                             distance=dist)
    
    return hp_tilde, hc_tilde

def a_sep(m1, m2, f):
    M = m1+m2
    a = (G*M*smtokg/(np.pi*f)**2)**(1/3)
    return a

def t_coal(m1, m2, f):
    M = m1+m2
    M_r = (m1*m2)/(m1+m2)
    t = 5/256*(c**5/(G**3*M**2*M_r*smtokg**3))*a_sep(m1,m2,f)**4
    return t

def t_bet_f(m1,m2,f1,f2):
    t1 = t_coal(m1, m2, f1)
    t2 = t_coal(m1, m2, f2)
    
    return t1-t2
    
def unite(a):
   # af = [complex(a[0], a[1])]
    af = []
    for i in range(len(a)-1):
        if i%2 == 0:
            af.append(complex(a[i],a[i+1]))
    af = np.array(af)
    return af

def divide(a):
    ad = []
    for i in range(len(a)):
        ad.append(a[i].real)
        ad.append(a[i].imag)
    ad = np.array(ad)
    return ad
#%
obs = pd.read_csv('/home/paolo/Dropbox/PhD/Python/wpd_datasets.csv')

def find_zeros(fun):
    l = len(fun)
    zeros = []
    for i in range(l-1):
        if fun[i+1]*fun[i] < 0:
            zeros.append(i)
    
    return zeros
#%
def f_from_wf(h_plus, h_cross):
    '''
    # phase = np.unwrap(np.arctan2(h_cross.data, h_plus.data))#.astype(real_same_precision_as(h_plus))

    phase = TimeSeries(phase, delta_t=h_plus.delta_t, epoch=h_plus.start_time)#), copy=False)
    # phase = utils.phase_from_polarizations(h_plus, h_cross)
    # phase = correct_phase(phase)
    # phase = phase_corrected
    
    freq = np.diff(phase) / ( 2 * np.pi * phase.delta_t )
    '''
    zeros = find_zeros(h_plus.data)
    f = []
    t_f = []
    # perc=1/50
    # '''
    times = np.array(h_plus.sample_times).round(8)
    for i in range(len(zeros)-1):
        # if i>perc*len(zeros):
            # print('completed %.1f %%'%(i/(len(zeros)-1)*100))
            # perc += 1/50
        # a = (times[zeros[i+1]]+times[zeros[i+1]+1])/2
        # b = 
        delta_t = times[zeros[i+1]]-times[zeros[i]]
        f.append(1/(2*delta_t))
        t_f.append(times[zeros[i]]+delta_t/2)
        # '''
    # start_time = phase.start_time + phase.delta_t / 2
    # f = TimeSeries(freq.astype(real_same_precision_as(h_plus)),
        # delta_t=hp.delta_t, epoch=start_time)
    
    return zeros, f, t_f

def f_tail(f2, f):
    f1 = 0.7*f2
    ta = (np.exp(((f2-f1)/(f-f1))+((f2-f1)/(f-f2)))+1)**-1
    
    return ta

def make_hp_tail(f2, h_tilde, factor=0.7):
    df = h_tilde.delta_f
    f1 = factor*f2
    fs = np.arange(f1+df, f2, df)
    ind_f2 = int(f2/df)
    ind_f1 = int(f1/df)
    
    f_exp = np.array([f_tail(f2, f) for f in fs])
    
    tail = np.concatenate((np.zeros(ind_f1), f_exp))
    l = int(len(h_tilde) - len(tail))
    tail = np.concatenate((tail, np.ones(l)))
    
    hp_t = h_tilde * tail
    
    return hp_t

def resize_h_tilde(h_tilde, fr, l_new, dt_new):
    '''
    l_p = 1/(dt*f_f-1/2)
    l_d = 1/(2*dt*f_f-1)
    
    if l_p % 2 == 0:
        l_new = l_p
    elif l_d % 2 == 1:
        l_new = l_d
    '''
    fr = np.array(fr)
    h_tilde = np.array(h_tilde)
    
    if l_new % 2 == 0:
        fr_new = np.arange((l_new/2)+1)/(l_new*dt_new)
    elif l_new % 2 == 1:
        fr_new = np.arange((l_new+1)/2)/(l_new*dt_new)

    if fr_new[-1] > fr[-1]:
        fr = np.append(fr, fr_new[-1])
        h_tilde = np.append(h_tilde, complex(0,0))
    
    h_f = interp1d(fr, h_tilde)
    
    h_tilde_new = h_f(fr_new)
    
    return h_tilde_new, fr_new

def AF_pm(f, M_lens, zL, y, MSD_fact=1):
    Ome = 8*np.pi*G*M_lens*(1+zL)*smtokg/c**3*f
    
    exp = cmath.exp((np.pi*Ome/4)+((1j*Ome/2)*cmath.log(Ome/2, np.e))) 
    g = gamma(1-(1j*Ome/2))
    hy = hyp1f1(1j*Ome/2,1,(1j*Ome*y**2)/2)
    
    af = exp*g*hy
    
    return complex(float(af.real), float(af.imag))#, Ome#, exp, g, hy

def theta_E(M, D_ls, D_d, D_s):
    er = np.sqrt(4*G*M*smtokg/(c**2)*D_ls/(D_d*D_s))
    
    return er

def AF_pm_go(f, M_lens, zL, y, lam=1, imag=0):
    # y = lam * y
    Ome = 8*np.pi*G*M_lens*(1+zL)*smtokg/c**3*f
    x_p, x_m = 1/2*(y+np.sqrt(y**2+4)), 1/2*(y-np.sqrt(y**2+4))
    mu_p, mu_m = ((1-(1/x_p)**4)**(-1)), ((1-(1/x_m)**4)**(-1))
    td_p = ((2*y**2+4-y*(y**2+4)**(1/2))/4 - np.log(abs((y+(y**2+4)**(1/2))/2)))
    td_m = ((2*y**2+4+y*(y**2+4)**(1/2))/4 - np.log(abs((y-(y**2+4)**(1/2))/2)))
    # et = np.exp(-1j*Ome*lam*(1-lam)/2*(1+zL)/(c**3)*4*G*M_lens*smtokg*y**2)
    M_red = (1+zL)*4*M_lens*smtokg/(c**3)
    # et = np.exp(-1j*Ome*lam*(1-lam)/2*G*M_red*(y**2))
    et = np.exp(-1j*Ome*lam*(1-lam)/2*(y**2))
    # af_1 = np.sqrt(abs(mu_p))*np.exp((-1)*1j*Ome*td_p)
    # af_2 = np.sqrt(abs(mu_m))*np.exp((-1)*((1j*Ome*td_m)-(1j*np.sign(Ome)*np.pi/2)))
    
    af_1 = 1/lam*np.sqrt(abs(mu_p))*np.exp(1j*Ome*lam*td_p)#*et
    af_2 = 1/lam*np.sqrt(abs(mu_m))*np.exp(((1j*Ome*lam*td_m)-(1j*np.sign(Ome)*np.pi/2)))#*et
    if imag == 'I':
        return af_1*np.exp(-1j*Ome*lam*td_p)
    if imag == "II":
        return af_2*np.exp(-1j*Ome*lam*td_m)
    else:
        return af_1 + af_2
   
def AF_PM(f, M_lens, zL, y, lam=1):
    Ome = 8*np.pi*G*M_lens*(1+zL)*smtokg/c**3*f
    q = (1j*Ome*lam)/2
    p = 2**(-1-q)*np.exp(q*lam*y**2)*1/lam*(-2*q)**(1+q)
    g = gamma(-q)
    # n = -1 + q
    # xx = -1/2*1j*Ome*lam*y**2
    # l = np.polynomial.laguerre.lagval(xx, np.eye(1, n+1, n)[0])
    # l = lpn(n, xx)
    h = hyp1f1(-(-1+q), 1, -q*y**2, maxterms=10**6)
    h = complex(h.real, h.imag)
    
    return p*g*h
   
def gps_datetime(gpstime):
    utc = datetime(1980, 1, 6) + timedelta(seconds=gpstime - (35 - 19))
    
    return utc

def make_amps(case, M_l, zL, y, hp_tilde, final_freq_index, dt=5*10**-5, ttype='wo', save=False, MSD_fact=1):
    file = 'amps_M'+str(int(M_l))+'_zl'+str(zL).replace('.','')+'_y'+str(round(y,2)).replace('.','')
    print('\n%s\n'%(file))
    if ttype == 'wo':
        AF = AF_pm
    elif ttype == 'go':
        AF = AF_pm_go
        
    l_t =len(hp_tilde) 
    l = (l_t-1)*2
    fr = np.arange(l/2+1)/(l*dt)
    amps = []
    pr=1/20
    ind_fl = final_freq_index
    for i in fr[:ind_fl]: 
        if i > fr[int((len(fr[:ind_fl]))*pr)]:
            # print('\r%.2f %%'%(pr*100), end='')
            print('\r'+str(round(pr*100, 2))+' % |'+'-'*int(pr*50)+' '*int(50*(1-pr))+'|', end='')
            pr = pr+(1/20)

        amps.append(AF(i, M_l, zL, y, MSD_fact))

    # print('\r%.2f %%'%(100), end='')
    print('\r100 % |'+'-'*50+'|', end='')    
    
    amps = np.array(amps)
    amps[0] = 0
    amps = np.concatenate((amps, np.zeros(l_t-len(amps))))
    
    if save == True:
        amps.tofile('/home/paolo/Desktop/waveform/case_'+case+'/PM/amps_new/'+file)
        
    return amps

def make_transform_with_C_c2r(array, fromPY=False):
    if fromPY:
        l = 2*len(array)
    else :
        l = 2*len(array)-2
    path = '/home/paolo/Documents/FT/'
    file_name = 'array'
    file_out = 'FT_array'
    
    array = divide(array)
    array.tofile(path+file_name)
    
    f = open(path+'file_c2r.c', 'w')
    f.write('#include <stdio.h>\n')
    f.write('#include <stdlib.h>\n')
    f.write('#include <complex.h>\n')
    f.write('#include <fftw3.h>\n')
    f.write('int main () {\n')
    f.write('        int N = '+str(l)+';\n')
    f.write('        fftw_complex *buffer = fftw_malloc(sizeof(double)*N+2);\n')
    f.write('        FILE *ptr, *fo;\n')
    f.write('        ptr = fopen("'+path+file_name+'", "rb");\n')
    f.write('        fread(buffer, sizeof(double), N+2, ptr);\n')
    f.write('        fftw_plan p;\n')
    f.write('        double *out = fftw_malloc(sizeof(double) * N);\n')
    f.write('        p = fftw_plan_dft_c2r_1d(N, buffer, out, FFTW_ESTIMATE);\n')
    f.write('        fftw_execute(p);\n')
    f.write('        fo = fopen("'+path+file_out+'", "w");\n')
    f.write('        fwrite(out, sizeof(double), N, fo);\n')
    f.write('        fclose(fo);\n')
    f.write('        fclose(ptr);\n')
    f.write('        fftw_destroy_plan(p);\n')
    f.write('        fftw_free(buffer); fftw_free(out);\n')
    f.write('}')
    
    f.close()
    
    os.system('gcc -o '+path+'o_ft_c2r '+path+'file_c2r.c -lfftw3 -lm')
    os.system(path+'./o_ft_c2r')
    
    output = np.fromfile(path+file_out)
    output /= len(output)
    
    return output

def make_transform_with_C_r2c(array):
    l = len(array)#2*len(array)-2
    path = '/home/paolo/Documents/FT/'
    file_name = 'td_array'
    file_out = 'fd_array'
    
    # array = divide(array)
    array.tofile(path+file_name)
    
    f = open(path+'file_r2c.c', 'w')
    f.write('#include <stdio.h>\n')
    f.write('#include <stdlib.h>\n')
    f.write('#include <complex.h>\n')
    f.write('#include <fftw3.h>\n')
    f.write('int main () {\n')
    f.write('        int N = '+str(l)+';\n')
    f.write('        double *buffer = fftw_malloc(sizeof(double)*N);\n')
    f.write('        FILE *ptr, *fo;\n')
    f.write('        ptr = fopen("'+path+file_name+'", "rb");\n')
    f.write('        fread(buffer, sizeof(double), N, ptr);\n')
    f.write('        fftw_complex *out;\n')
    f.write('        fftw_plan p;\n')
    f.write('        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N+2);\n')
    f.write('        p = fftw_plan_dft_r2c_1d(N, buffer, out, FFTW_ESTIMATE);\n')
    f.write('        fftw_execute(p);\n')
    f.write('        fo = fopen("'+path+file_out+'", "w");\n')
    f.write('        fwrite(out, sizeof(double), N+2, fo);\n')
    f.write('        fclose(fo);\n')
    f.write('        fclose(ptr);\n')
    f.write('        fftw_destroy_plan(p);\n')
    f.write('        fftw_free(buffer); fftw_free(out);\n')
    f.write('}')
    
    f.close()
    
    os.system('gcc -o '+path+'o_ft_r2c '+path+'file_r2c.c -lfftw3 -lm')
    os.system(path+'./o_ft_r2c')
    
    output = np.fromfile(path+file_out)
    output = unite(output)
    
    return output


def plot_l(label, tL, hp_L, align = 'no', pu=0, pl=0, n_roll = 0, color = 'orangered', ls = None):
    
    
    if align == 'si' :
        # pu = np.where(hp==max(hp))[0][0]
        # pl = np.where(hp_L==max(hp_L))[0][0]
        plt.plot(tL, np.roll(hp_L, pu-pl), linestyle=ls, c=color, label=label)
    elif align == 'no' :
        plt.plot(tL, hp_L, c=color, label=label)
    elif align == 'custom' :
        plt.plot(tL, np.roll(hp_L, int(n_roll)), c=color, label=label)
        
def make_hpL(amps_file, unlensed_fd_file, fd=False, fromPY=False):
    if type(amps_file) == np.str :
        amps = np.fromfile(amps_file, np.complex64)
    elif type(amps_file) == np.ndarray : 
        amps = amps_file
    '''
    if len(amps) == 2157:
        hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74'
        hp_t_fn = '/home/paolo/Desktop/waveform/case_h/hp_fd_H74'
        l = 2157032
    elif len(amps) == 1262:
        hp_fn = '/home/paolo/Desktop/waveform/case_h/hp_td_H74_mt'
        hp_t_fn = '/home/paolo/Desktop/waveform/case_h/hp_fd_H74_mt'
        l = 3157032
    '''
    if type(unlensed_fd_file) == np.str :
        hp_t = np.fromfile(unlensed_fd_file)
        hp_t = unite(hp_t)
    elif type(unlensed_fd_file) == np.ndarray:
        hp_t = unlensed_fd_file
        
    amps = np.concatenate((amps, np.zeros(len(hp_t)-len(amps))))
    hp_L_t = hp_t * np.conj(amps)
    
    if not fd:
        hp_L = make_transform_with_C_c2r(hp_L_t, fromPY)
        # hp_L /= len(hp_L)
        return hp_L
    else:
        return hp_L_t
    
def inner_p(a, fr_a, b, fr_b, sn, fr_sn, fr_min, fr_max):
    fa_r = interp1d(fr_a, np.real(a))
    fb_r = interp1d(fr_b, np.real(b))

    fa_i = interp1d(fr_a, np.imag(a))
    fb_i = interp1d(fr_b, np.imag(b))

    fs = interp1d(fr_sn, sn)

    def fa(x) : return fa_r(x)+1j*fa_i(x)
    def fb(x) : return fb_r(x)+1j*fb_i(x)
    def fp(x): return (fa(x)*np.conj(fb(x)))/fs(x)
    
    ll = fr_min#max(fr_a[0], fr_b[0], fr_sn[0])
    ul = fr_max #275#155#min(fr_a[-1], fr_b[-1], fr_sn[fr_sn.index[-1]])
    '''
    res = 4*quadpy.quad(fp, ll, ul, limit=1000)
    return np.real(res[0])
    '''
    df = fr_a[1]-fr_a[0]
    fr = np.arange(ll, ul, df)
    in_i = np.where(fr>ll)[0][0]
    in_f = np.where(fr<ul)[0][-1]
    res = np.sum(fp(fr)*df)
    # res = quadpy.quad(fp, ll, ul, limit=10000)[0]
    # res_1 = (a[in_i:in_f]*np.conj(b[in_i:in_f]))/fs(fr[in_i:in_f])
    # res = 4*np.sum(res_1*df)
    return 4*res

# from Hezquiaga 

def Fp(theta,phi,Psi): #(57) of https://arxiv.org/pdf/0903.0338.pdf
    return 0.5*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.cos(2*Psi) - np.cos(theta)*np.sin(2*phi)*np.sin(2*Psi)

def Fx(theta,phi,Psi): #(57) of https://arxiv.org/pdf/0903.0338.pdf
    return 0.5*(1 + np.cos(theta)**2)*np.cos(2*phi)*np.sin(2*Psi) + np.cos(theta)*np.sin(2*phi)*np.cos(2*Psi)

def mchirp(m1,m2):
    return np.power(m1*m2,3./5.)/np.power(m1+m2,1./5.)


def f_PhenomA(indx,m1,m2):
    a = np.array([2.9740e-1, 5.9411e-1, 5.0801e-1, 8.4845e-1])
    b = np.array([4.4810e-2, 8.9794e-2, 7.7515e-2, 1.2848e-1])
    c = np.array([9.5560e-2, 1.9111e-1, 2.2369e-2, 2.7299e-1])
    tM = (m1 + m2)*TSUN 
    eta = m1*m2/((m1+m2)**2.)
    num = a[indx]*eta**2 + b[indx]*eta + c[indx]
    return num/(np.pi*tM)

def f_fin(fi,m1,m2, T=10):#,T,baseline): 
    #fi in Hz
    #M in Msun in *detector* frame
    #T of mission in yr
    ffin = f_PhenomA(3,m1,m2)
    '''
    if T < t_merge(fi,Mc)/YEAR:
        tM = Mc*TSUN 
        ffin = (1./(8.*np.pi*tM))*(-(T*YEAR)/(5.*tM) + (1./(8.*np.pi*tM*fi))**(8./3.))**(-3./8.)
    else:
        ffin = f_PhenomA(3,m1,m2)
        if baseline == 'space':
            ffin = min(ffin,2.)
    '''
    return ffin # Hz
vf_fin = np.vectorize(f_fin)
