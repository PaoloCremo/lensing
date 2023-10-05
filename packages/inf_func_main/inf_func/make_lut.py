# packages
import os
path_home = os.path.expanduser('~')
import sys

import bilby
import numpy as np
from scipy.special import gamma
from mpmath import hyp1f1, re, im
from scipy.interpolate import interpn
import datetime

from inf_func import lut

# constants

c = 299_792_459       #m/s
G = 6.67408*10**(-11) #m^3 kg^-1 s^-2
smtokg=1.98847*10**30 #kg/M_\odot
TSUN    = 4.92549232189886339689643862e-6 # mass of sun in seconds (G=C=1)i

# define functions

hyv = np.vectorize(hyp1f1)
rev, imv = np.vectorize(re), np.vectorize(im)
compv = np.vectorize(complex)

def AF_pm_GO(f, Mlr, y, lam=1, imag=0):
    
    Ome = 8*np.pi*G*Mlr*smtokg/c**3*f
    x_p, x_m = 1/2*(y+np.sqrt(y**2+4)), 1/2*(y-np.sqrt(y**2+4))
    mu_p, mu_m = ((1-(1/x_p)**4)**(-1)), ((1-(1/x_m)**4)**(-1))
    td_p = ((2*y**2+4-y*(y**2+4)**(1/2))/4 - np.log(abs((y+(y**2+4)**(1/2))/2)))
    td_m = ((2*y**2+4+y*(y**2+4)**(1/2))/4 - np.log(abs((y-(y**2+4)**(1/2))/2)))
    M_red = 4*Mlr*smtokg/(c**3)
    et = np.exp(-1j*Ome*lam*(1-lam)/2*(y**2))
    af_1 = 1/lam*np.sqrt(abs(mu_p))*np.exp(1j*Ome*lam*td_p)#*et
    af_2 = 1/lam*np.sqrt(abs(mu_m))*np.exp(((1j*Ome*lam*td_m)-(1j*np.sign(Ome)*np.pi/2)))#*et
    if imag == 'I':
        return af_1*np.exp(-1j*Ome*lam*td_p)
    if imag == "II":
        return af_2*np.exp(-1j*Ome*lam*td_m)
    else:
        return af_1 + af_2

def AF_PM2(f, M_red, y, lam=1):
    '''
    function to calculate the amplification factor
    f:     frequency
    M_red: redshifted (observer frame) mass of the lens
    y:     impact parameter (dimensionless)
    lam:   MSD parameter - lam=1 no transformation
    '''
    Ome = 8*np.pi*G*M_red*smtokg/c**3*f
    q = (1j*Ome*lam)/2
    p = 2**(-1-q)*np.exp(q*lam*y**2)*1/lam*(-2*q)**(1+q)
    g = gamma(-q)
    h = hyv(-(-1+q), 1, -q*y**2, maxterms=10**6)
    h = compv(rev(h), imv(h))

    return p*g*h

def AF_int(matrix, f, Ml_red, y, lam):
    '''
    function to interpolate the LUT with parameters from bilby
    '''
    return interpn((lut.f_a, lut.Mlr_a, lut.y_a, lut.lam_a), matrix, (f, Ml_red, y, lam))

def AF_int_bilby_pipe(f, Ml_red, y, lam):
    '''
    function to interpolate the LUT with parameters from bilby
    '''

    return interpn((lut.f_a, lut.Mlr_a, lut.y_a, lut.lam_a), lut.matrix, (f, Ml_red, y, lam))



def my_lal_binary_black_hole(frequency_array, mass_1, mass_2, luminosity_distance, 
                             a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, 
                             phase, Ml_r, y, lam, **kwargs):
    '''
    function to create the waveform to feed the MCMC
    '''
    waveform_kwargs = dict(waveform_approximant='IMRPhenomXP', reference_frequency=50.0, 
                           minimum_frequency=20.0, maximum_frequency=frequency_array[-1], 
                           catch_waveform_errors=False, pn_spin_order=-1, pn_tidal_order=-1, 
                           pn_phase_order=-1, pn_amplitude_order=0)
    waveform_kwargs.update(kwargs)
    wf = bilby.gw.source._base_lal_cbc_fd_waveform( frequency_array=frequency_array, 
                                                   mass_1=mass_1, mass_2=mass_2,
                                                   luminosity_distance=luminosity_distance, 
                                                   theta_jn=theta_jn, phase=phase, a_1=a_1, 
                                                   a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, 
                                                   phi_12=phi_12, phi_jl=phi_jl, 
                                                   **waveform_kwargs)
    
    lf = np.where(frequency_array.round(0)==19)[0][-1]           # find index for f_min
    f_max = round(f_fin(20., mass_1, mass_2), 0)                 # find f_max
    uf = np.where(frequency_array.round(1)==f_max)[0][0]         # find index for f_max
    AF = AF_int(matrix, frequency_array[lf:uf], Ml_r, y, lam)    # compute AF in [f_min:f_max] with the LUT
    AF = np.concatenate((np.zeros(lf), AF))                      # add zeros for f<f_min
    AF = np.concatenate((AF, np.zeros(len(frequency_array)-uf))) # add zeros for f>f_max

    return dict(plus=wf['plus']*np.conj(AF), cross=wf['cross']*np.conj(AF))

def my_lal_binary_black_hole_bilby_pipe(frequency_array, mass_1, mass_2, luminosity_distance,
                             a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
                             phase, Ml_r, y, lam, **kwargs):
    '''
    function to create the waveform to feed the MCMC
    '''
    waveform_kwargs = dict(waveform_approximant='IMRPhenomXP', reference_frequency=50.0,
                           minimum_frequency=20.0, maximum_frequency=frequency_array[-1],
                           catch_waveform_errors=False, pn_spin_order=-1, pn_tidal_order=-1,
                           pn_phase_order=-1, pn_amplitude_order=0)
    waveform_kwargs.update(kwargs)
    wf = bilby.gw.source._base_lal_cbc_fd_waveform( frequency_array=frequency_array,
                                                   mass_1=mass_1, mass_2=mass_2,
                                                   luminosity_distance=luminosity_distance,
                                                   theta_jn=theta_jn, phase=phase, a_1=a_1,
                                                   a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2,
                                                   phi_12=phi_12, phi_jl=phi_jl,
                                                   **waveform_kwargs)

    lf = np.where(frequency_array.round(0)==19)[0][-1]           # find index for f_min
    f_max = round(f_fin(20., mass_1, mass_2), 0)                 # find f_max
    uf = np.where(frequency_array.round(1)==f_max)[0][0]         # find index for f_max
    AF = AF_int_bilby_pipe(frequency_array[lf:uf], Ml_r, y, lam)    # compute AF in [f_min:f_max] with the LUT
    AF = np.concatenate((np.zeros(lf), AF))                      # add zeros for f<f_min
    AF = np.concatenate((AF, np.zeros(len(frequency_array)-uf))) # add zeros for f>f_max

    return dict(plus=wf['plus']*np.conj(AF), cross=wf['cross']*np.conj(AF))


def my_lal_binary_black_hole2(frequency_array, mass_1, mass_2, luminosity_distance, 
                             a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, 
                             phase, Ml_r, y, lam, **kwargs):
    '''
    function to create the waveform to feed the MCMC
    '''
    waveform_kwargs = dict(waveform_approximant='IMRPhenomXP', reference_frequency=50.0, 
                           minimum_frequency=20.0, maximum_frequency=frequency_array[-1], 
                           catch_waveform_errors=False, pn_spin_order=-1, pn_tidal_order=-1, 
                           pn_phase_order=-1, pn_amplitude_order=0)
    waveform_kwargs.update(kwargs)
    wf = bilby.gw.source._base_lal_cbc_fd_waveform( frequency_array=frequency_array, 
                                                   mass_1=mass_1, mass_2=mass_2,
                                                   luminosity_distance=luminosity_distance, 
                                                   theta_jn=theta_jn, phase=phase, a_1=a_1, 
                                                   a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, 
                                                   phi_12=phi_12, phi_jl=phi_jl, 
                                                   **waveform_kwargs)
    
    lf = np.where(frequency_array.round(0)==19)[0][-1]           # find index for f_min
    f_max = round(f_fin(20., mass_1, mass_2), 0)                 # find f_max
    uf = np.where(frequency_array.round(1)==f_max)[0][0]         # find index for f_max
    AF = AF_PM2(frequency_array[lf:uf], Ml_r, y, lam)            # compute AF in [f_min:f_max] 
    AF = np.concatenate((np.zeros(lf), AF))                      # add zeros for f<f_min
    AF = np.concatenate((AF, np.zeros(len(frequency_array)-uf))) # add zeros for f>f_max

    return dict(plus=wf['plus']*np.conj(AF), cross=wf['cross']*np.conj(AF))

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

def make_table(f_array, Ml_array, y_array, lam_array, file_path = path_home+'/matrix'):
    # zL = zL_array
    m_re = np.zeros((len(f_array), len(Ml_array), len(y_array), len(lam_array)))
    m_im = np.zeros((len(f_array), len(Ml_array), len(y_array), len(lam_array)))
    # m_re = np.zeros((len(Ml_array), len(y_array), len(lam_array), len(f_array)))
    # m_im = np.zeros((len(Ml_array), len(y_array), len(lam_array), len(f_array)))
    
    comp = np.vectorize(complex)
    
    ll = len(lam_array)
    ly = len(y_array)
    lm = len(Ml_array)
    lp1 = ly*lm
    lp2 = lp1*ll
    AFv = np.vectorize(AF_PM2)
    
    now = datetime.datetime.now()
    for a,lam in enumerate(lam_array):
        for b,y in enumerate(y_array):
            for c,Ml in enumerate(Ml_array):
                '''
                m_re[c][b][a] = np.concatenate(([0], np.real(AFv(f_array[1:], Ml, y, lam))))
                m_im[c][b][a] = np.concatenate(([0], np.imag(AFv(f_array[1:], Ml, y, lam))))
                # m_re[1:][c][b][a] = np.real(AFv(f_array[1:], Ml, y, lam))
                # m_im[1:][c][b][a] = np.imag(AFv(f_array[1:], Ml, y, lam))
                '''
                m_re[0][c][b][a] = 0
                m_im[0][c][b][a] = 0
                for d,f in enumerate(f_array[1:]):
                    m_re[d+1][c][b][a] = np.real(AF_PM2(f, Ml, y, lam))
                    m_im[d+1][c][b][a] = np.imag(AF_PM2(f, Ml, y, lam))
            
            with open('output_m.txt', 'w') as f:
                pct = ((a)*lp1+(b)*lm+(c+1))/lp2*100
                td = datetime.datetime.now() - now
                total = td*100/pct
                rem = (total-td)*1.5
                print('lambda: %i of %i | y: made  %i of %i\ncompleted: %.3f %%\n'%(a, len(lam_array), b+1, len(y_array), pct),  file=f)
                print('elapsed: {:>17s}\ntotal exp: {:>15s}\nremaining: {:>15s}'.format(str(td), str(rem+td), str(rem)), file=f) 

            if b%10==0:
                m = comp(m_re, m_im)
                m.tofile(file_path)
    
    m = comp(m_re, m_im)
    m.tofile(file_path)

    return None

noe = 0
def run_af(f, Ml, y, lam, dm):
    global noe
    while True:
        try:
            af = AF_PM2(f, Ml, y, lam)
        except OverflowError:
            if noe == 0:
                with open('exceptions.txt', 'a') as ff:
                    print('ex_n: f; Ml; y; lam\n', file=ff)
            noe += 1
            with open('exceptions.txt', 'a') as ff:
                print('{}: {}; {}; {:.2f}; {:.2f}'.format(noe, f, Ml, y, lam), file=ff)
            Ml -= dm
            continue
        break

    return af

def run_af_mix(f, Ml, y, lam, dm):
    global noe
    while True:
        try:
            af = AF_PM2(f, Ml, y, lam)
        except OverflowError:
            if noe == 0:
                with open('exceptions.txt', 'a') as ff:
                    print('ex_n: f; Ml; y; lam\n', file=ff)
            noe += 1
            with open('exceptions.txt', 'a') as ff:
                print('{}: {}; {}; {:.2f}; {:.2f}'.format(noe, f, Ml, y, lam), file=ff)
            # af = AF_pm_GO(f, Ml, y, lam)
            af = 1. 
            continue
        break

    return af


def make_table2(f_array, Ml_array, y_array, lam_array, file_path = path_home+'/matrix', outfile='output_m.txt'):
    # zL = zL_array
    m_re = np.zeros((len(f_array), len(Ml_array), len(y_array), len(lam_array)))
    m_im = np.zeros((len(f_array), len(Ml_array), len(y_array), len(lam_array)))
    # m_re = np.zeros((len(Ml_array), len(y_array), len(lam_array), len(f_array)))
    # m_im = np.zeros((len(Ml_array), len(y_array), len(lam_array), len(f_array)))
    
    comp = np.vectorize(complex)
    
    ll = len(lam_array)
    ly = len(y_array)
    lm = len(Ml_array)
    lp1 = ly*lm
    lp2 = lp1*ll
    AFv = np.vectorize(AF_PM2)
    
    now = datetime.datetime.now()
    dm = Ml_array[-1]-Ml_array[-2] 
    for a,lam in enumerate(lam_array):
        for b,y in enumerate(y_array):
            for c,Ml in enumerate(Ml_array):
                '''
                m_re[c][b][a] = np.concatenate(([0], np.real(AFv(f_array[1:], Ml, y, lam))))
                m_im[c][b][a] = np.concatenate(([0], np.imag(AFv(f_array[1:], Ml, y, lam))))
                # m_re[1:][c][b][a] = np.real(AFv(f_array[1:], Ml, y, lam))
                # m_im[1:][c][b][a] = np.imag(AFv(f_array[1:], Ml, y, lam))
                '''
                m_re[0][c][b][a] = 0
                m_im[0][c][b][a] = 0
                for d,f in enumerate(f_array[1:]):
                    af_v = run_af_mix(f, Ml, y, lam, dm)
                    m_re[d+1][c][b][a] = np.real(af_v)
                    m_im[d+1][c][b][a] = np.imag(af_v)
            
            with open(outfile, 'w') as ff:
                pct = ((a)*lp1+(b)*lm+(c+1))/lp2*100
                td = datetime.datetime.now() - now
                total = td*100/pct
                rem = (total-td)*1.5
                print('lambda: %i of %i | y: made  %i of %i\ncompleted: %.3f %%\n'%(a, len(lam_array), b+1, len(y_array), pct),  file=ff)
                print('elapsed: {:>17s}\ntotal exp: {:>15s}\nremaining: {:>15s}'.format(str(td), str(rem+td), str(rem)), file=ff) 

            if b%10==0:
                m = comp(m_re, m_im)
                m.tofile(file_path)
    
    m = comp(m_re, m_im)
    m.tofile(file_path)

    return None


def main():
    make_table2(lut.f_a6, lut.Mlr_a6, lut.y_a6, lut.lam_a6)