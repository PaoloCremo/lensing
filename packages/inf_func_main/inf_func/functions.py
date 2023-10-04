# packages
import bilby
import numpy as np
from scipy.interpolate import interpn

from inf_func import lut

# constants

c = 299_792_459       #m/s
G = 6.67408*10**(-11) #m^3 kg^-1 s^-2
smtokg=1.98847*10**30 #kg/M_\odot
TSUN    = 4.92549232189886339689643862e-6 # mass of sun in seconds (G=C=1)i

# define functions

def AF_int_bilby_pipe(f, Ml_red, y, lam, matrix_n):
    '''
    function to interpolate the LUT with parameters from bilby

    matrix_n : LUT for computing amplification factor. 
               Available n = [1, 3, 4, 5]
    '''
    available = [1,3,4,5]
    if matrix_n not in available:
        raise TypeError('"matrix_n" must be one of {}'.format(available))

    elif matrix_n == 1:
        return interpn((lut.f_a, lut.Mlr_a, lut.y_a, lut.lam_a), lut.matrix, (f, Ml_red, y, lam))
    elif matrix_n == 3: 
        return interpn((lut.f_a3, lut.Mlr_a3, lut.y_a3, lut.lam_a3), lut.matrix3, (f, Ml_red, y, lam))
    elif matrix_n == 4:
        return interpn((lut.f_a4, lut.Mlr_a4, lut.y_a4, lut.lam_a4), lut.matrix4, (f, Ml_red, y, lam))
    elif matrix_n == 5:
        return interpn((lut.f_a5, lut.Mlr_a5, lut.y_a5, lut.lam_a5), lut.matrix5, (f, Ml_red, y, lam))

def my_lal_binary_black_hole_bilby_pipe(frequency_array, mass_1, mass_2, luminosity_distance,
                             a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
                             phase, Ml_r, y, lam, matrix_n, **kwargs):
    '''
    function to create the waveform to feed the MCMC
    '''
    waveform_kwargs = dict(waveform_approximant='IMRPhenomXPHM', reference_frequency=50.0,
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
    AF = AF_int_bilby_pipe(frequency_array[lf:uf], Ml_r, y, lam, matrix_n) # compute AF in [f_min:f_max] with the LUT
    AF = np.concatenate((np.zeros(lf), AF))                      # add zeros for f<f_min
    AF = np.concatenate((AF, np.zeros(len(frequency_array)-uf))) # add zeros for f>f_max

    return dict(plus=wf['plus']*np.conj(AF), cross=wf['cross']*np.conj(AF))

def my_lal_binary_black_hole_bilby_pipe_m4(frequency_array, mass_1, mass_2, luminosity_distance,
                             a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
                             phase, Ml_r, y, lam, **kwargs):
    return my_lal_binary_black_hole_bilby_pipe(frequency_array, mass_1, mass_2, luminosity_distance,
                             a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
                             phase, Ml_r, y, lam, matrix_n=4, **kwargs)

def my_lal_binary_black_hole_bilby_pipe_m5(frequency_array, mass_1, mass_2, luminosity_distance,
                             a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
                             phase, Ml_r, y, lam, **kwargs):
    return my_lal_binary_black_hole_bilby_pipe(frequency_array, mass_1, mass_2, luminosity_distance,
                             a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn,
                             phase, Ml_r, y, lam, matrix_n=5, **kwargs)
    


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
