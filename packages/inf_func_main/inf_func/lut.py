import os
path_home = os.path.expanduser('~')

import numpy as np

# frequencies
f1 = np.arange(0, 100, 1)
f2 = np.arange(100, 1000, 10)
f3 = np.arange(1000, 5000, 50)
# arrays
f_a = np.concatenate((f1, f2, f3)) # frequencies
Mlr_a = np.arange(88, 132.1, 4)    # observer lens mass
y_a = np.arange(0.5, 2.01, 0.04)   # impact parameter
lam_a = np.arange(0.8, 1.21, 0.03) # MSD parameter

'''
matrix = np.fromfile(path_home+'/myscratch/matrices/matrix', dtype=np.complex)   # read LUT
matrix = np.reshape(matrix, (len(f_a), len(Mlr_a), len(y_a), len(lam_a)))  # reshape LUT to arrays

# matrix 3
# frequencies
f1_3 = np.arange(0, 100, 1)
f2_3 = np.arange(100, 1000, 10)
f3_3 = np.arange(1000, 5000, 50)
# arrays
f_a3 =   np.concatenate((f1_3, f2_3, f3_3)) # frequencies
Mlr_a3 = np.arange(60, 161.1, 4)            # observer lens mass
y_a3 =   np.arange(0.06, 1.94, 0.04)        # impact parameter
lam_a3 = np.arange(0.5, 1.53, 0.03)         # MSD parameter

matrix3 = np.fromfile(path_home+'/myscratch/matrices/matrix3', dtype=np.complex)   # read LUT
matrix3 = np.reshape(matrix3, (len(f_a3), len(Mlr_a3), len(y_a3), len(lam_a3)))  # reshape LUT to arrays

def conc_matrix(lf, lm, ly):
    path_folder = path_home+'/myscratch/matrices/'
    li = os.listdir(path_folder)
    li.sort()
    for n,file in enumerate(li):
        mp = np.fromfile(path_folder+file, dtype=complex)
        mp = np.reshape(mp, (lf, lm, ly, 1))
        if n==0:
            matrix = mp
        else:
            matrix = np.concatenate((matrix,mp), axis=3)
    return matrix


'''

# matrix4

#frequencies
f1_4 = np.arange(0, 100, 1)
f2_4 = np.arange(100, 1000, 10)
f3_4 = np.arange(1000, 5000, 25)

#masses
m1_4 = np.arange(50, 150, 4)
m2_4 = np.arange(150, 500, 10)
m3_4 = np.arange(500, 1001, 20)

# arrays
f_a4   = np.concatenate((f1_4, f2_4, f3_4))     # frequencies
Mlr_a4 = np.concatenate((m1_4, m2_4, m3_4))         # observer lens mass
y_a4   = np.arange(0.06, 2.54, 0.04)      # impact parameter
lam_a4 = np.arange(0.5, 1.53, 0.03)       # MSD parameter

matrix4 = np.fromfile(path_home+'/data/matrices/matrix4', dtype=complex)   # read LUT
matrix4 = np.reshape(matrix4, (len(f_a4), len(Mlr_a4), len(y_a4), len(lam_a4)))  # reshape LUT to arrays

# matrix 5

# frequencies
f1_5 = np.arange(0, 100, 1.5)
f2_5 = np.arange(100, 471, 2.5)
#masses
m1_5 = np.arange(200, 500, 10)
m2_5 = np.arange(500, 2000, 20)
m3_5 = np.arange(2000, 10001, 100)

# arrays
f_a5   = np.concatenate((f1_5, f2_5))         # frequencies
Mlr_a5 = np.concatenate((m1_5, m2_5, m3_5))     # observer lens mass
y_a5   = np.arange(0.06, 2.54, 0.04)      # impact parameter
lam_a5 = np.arange(0.5, 1.53, 0.03)       # MSD parameter

matrix5 = np.fromfile(path_home+'/data/matrices/matrix5', dtype=complex)   # read LUT
matrix5 = np.reshape(matrix5, (len(f_a5), len(Mlr_a5), len(y_a5), len(lam_a5)))  # reshape LUT to arrays
