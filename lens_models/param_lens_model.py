#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 18:03:31 2021

@author: paolo
"""

lmod = ['SIS', 'NIS', 'NFW', 'NFW_2', 'GNFW_0', 'GNFW_2', 'Hernquist', 'Burkert']
cols = ['red', 'yellow', 'k', 'k', 'green','green', 'b', 'cyan']
lss = [None, None, 'dashed', None, 'dashed', None, None, None]
ys = ['1', '1', '001061', '002804', '00287', '00634', '002677', '002802']
#%%
df = pd.DataFrame(data=[cols, lss, ys], columns = lmod)
df.to_csv('/home/paolo/Desktop/waveform/lens_model/codes/parameters.txt', index=False)
#%%
df = pd.read_csv('/home/paolo/Desktop/waveform/lens_model/codes/parameters.txt', index_col=False)

#%%
lmod = ['SIS', 'NFW', 'MoG', 'MoG', 'NFW_2', 'MoG_2', 'MoG_2']
cols = ['red', 'k', 'k', 'k', 'b', 'b', 'b']
lss = [None, None, 'dotted', 'dashed', None, 'dotted', 'dashed']
ys = ['1', '001061', '001061', '001061','002804','002804','002804']
#%%
df = pd.DataFrame(data=[cols, lss, ys], columns = lmod)
df.to_csv('/home/paolo/Desktop/waveform/lens_model/codes/parameters_MoG.txt', index=False)
#%%
df = pd.read_csv('/home/paolo/Desktop/waveform/lens_model/codes/parameters_MoG.txt', index_col=False)
