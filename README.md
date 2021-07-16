# files

### folder with models name (SIS, NFW, NFW-2 & GNFW_2)
Contain Mathematica script to calculate amplification factor

### packages
- `pyCBC_function.py` basic functions to calculate waveforms and lensing
- `LISA.py` LISA sensitivities from [here](https://github.com/eXtremeGravityInstitute/LISA_Sensitivity/tree/1)

### MSD
Codes used in the paper [arXiv:210407055](https://arxiv.org/abs/2104.07055)

- `MSD_2.py` general code to compute lensing 
- `MSD_matrix_2.py` code to plot matrix of SNRs
- `SNR_regimes.py` code to compute SNR calculation and plot
- `SNR_regimes_2.py` edited version of *"SNR_regimes.py"*

### microlensing
- `multi_lenses.py` code to compute lensing of sigle signal by multiple lenses

### lens_models
- `SIS_NFW.py`
- `SNR.py` calculate SNR
- `lens_models.py` 
- `lens_profiles.py` compute and plot lens mass profiles
- `param_lens_model.py` contain data about lens models
- `phase.py` compute phase calculations and plots
- `plots.py` compute different kind of plots

# lens models parameters

**amps files** </br>
*..._FP*    means for phase  and with the normalization factor phi_m </br>
*..._norm*  means long array and with the normalization factor phi_m </br>
*..._FP_nn* means for phase  and WITHOUT normalization </br>
*... * means long array and WITHOUT normalization </br>


m1 = m2 = 5*10^7</br> 
zS = 1</br>
dt = 1</br>
N = 2157032</br>
N_mt = 3157032</br>

### models parameters

#### M_L = 10^9 (rest frame)

zL = 0.5</br>
theta_E_SIS = 2.62394 * 10^-7
		
lens |sigma [km/s]| r @ y=1 : theta_E*DL [m]
-----|------------|-------------------------
SIS | 66.256 | 9.62299*10^18  
NIS | 66.10 | 9.62299*10^18 
 *lm* |**rho_s [kg/m^3]** | **rs [m]**
NFW | 7.954495462000312*10^-22 | 9.065935*10^20
NFW_2 | 7.954495*10^-22 | 3.4313*10^20
Hernquist | 2.87*10^-21 | 3.5948*10^20
GNFW 0 | 1.10*10^-21 | 3.3531*10^20
GNFW 2 | 4.938784307737533*10^-22 | 1.517823*10^20
BUR | 7.323121259691678*10^-22 | 3.43436*10^20
 
------------------

- **ySIS = 1.00000** (*red*) 
- yNIS = 1.00000    
- yNFW = 0.01061 (*blue*)  
- yNF2 = 0.02804 (* df_clr[9] , orangered (zs=0.55)  *)
- yHER = 0.02677    
- yGN0 = 0.02870    
- yGN2 = 0.06340 (*green*)
- yBUR = 0.02802

------------------

- **ySIS = 3.14100** (*red*)   
- yNFW = 0.02135 (*blue*)  
- yNF2 = 0.05640 (* df_clr[9] , orangered (zs=0.55)  *)   
- yGN2 = 0.12749 (*green*)

------------------

- **ySIS = 0.10000**
- yNIS = 0.10000
- yNFW = 0.00106
- yNF2 = 0.00280
- yHER = 0.00268
- yGN0 = 0.00287
- yGN2 = 0.00634
- yBUR = 0.00280

------------------

- **ySIS - 2.00000**
- yNIS - 2.00000
- yNFW - 0.02123
- yNF2 - 0.05609
- yHER - 0.05354
- yGN0 - 0.05740
- yGN2 - 0.12680
- yBUR - 0.05604

------------------

#### M_L = 5*10^8 (rest frame)

zL = 0.5</br>
theta_E_SIS = 1.86 * 10^-7
		
lens |sigma [km/s]| r @ y=1 : theta_E*DL [m]
-----|------------|-------------------------
SIS |         55.71 	|     	                6.80448*10^18  
NIS     |     55.58      |                     6.80448*10^18 
(*lm*) | **rho_s	[kg/m^3]**|		**rs [m]**
NFW |          7.954495462000312*10^-22   |     8.39304*10^20
NFW_2   |     7.954495*10^-22    |             2.88202*10^20
Hernquist |	2.87*10^-21 |                     3.32636*10^20
GNFW 0	| 1.10*10^-21     |                2.82031*10^20
GNFW 2 |	4.938784307737533*10^-22 |	1.2713*10^20
BUR |         7.323121259691678*10^-22 |       2.8819*10^20

------------------

- **ySIS - 1.00000**
- yNIS - 1.00000
- yNFW - 0.00811
- yNF2 - 0.02361
- yHER - 0.02046
- yGN0 - 0.02413
- yGN2 - 0.05352
- yBUR - 0.02356

------------------

- **ySIS = 0.10000**
- yNIS = 0.10000
- yNFW = 0.00081
- yNF2 = 0.00236
- yHER = 0.00205
- yGN0 = 0.00241
- yGN2 = 0.00535
- yBUR = 0.00236

------------------
