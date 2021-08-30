# -*- coding: utf-8 -*-
"""
Created on Tue May 18 19:54:51 2021

@author: pahi9557
"""

#When it actually comes down to detectibility in our pairs, what we're effectivley 
#measuring is 

import numpy as np
import matplotlib.pyplot as plt
# import pandas
import pandas as pd


class CTTS:
    def __init__(self,a,b,c,d,e,f,g,h,i,j,k):
        self.name = a
        self.accretion_rate = b
        self.num_of_5_sigma_expected_si4 = c
        self.num_of_5_sigma_expected_si4_err_up = d
        self.num_of_5_sigma_expected_si4_err_down = e
        self.num_of_3_sigma_expected_si4 = f
        self.num_of_3_sigma_expected_si4_err_up = g
        self.num_of_3_sigma_expected_si4_err_down = h
        self.AD_Leo_scale_ratio = i
        self.CTTS_RMS = j
        self.F_q_wtts = k
        
# TW_Hya = CTTS('TW Hya', 1.8e-9,4,0,0,4,1,2,0.15,0.044e-12,2.2e-14)
# GM_Aur = CTTS('GM Aur', 9.6e-9,10,4,3,9,9,4,0.15,0.019e-12,2.2e-14)
Recx_15 = CTTS('RECX 15', 8e-10, 32,16,10,64,12,21,0.11,0.0044e-12,1.6e-14)
# Recx_11 = CTTS('RECX 11', 1.7e-10, 24,16,6,24,24,14,0.11,0.0063e-12,1.6e-14)

#new flare numbers after hilton et al appraoch
TW_Hya = CTTS('TW Hya', 1.8e-9,4,0,0,1,1,1,0.15,0.044e-12,2.2e-14)
GM_Aur = CTTS('GM Aur', 9.6e-9,10,4,3,11,10,4,0.15,0.019e-12,2.2e-14)
Recx_11 = CTTS('RECX 11', 1.7e-10, 24,16,6,32,43,16,0.11,0.0063e-12,1.6e-14)






Accretion_rates = [TW_Hya.accretion_rate, GM_Aur.accretion_rate, Recx_11.accretion_rate]
names = [TW_Hya.name, GM_Aur.name, Recx_11.name]
scale_ratios = [TW_Hya.AD_Leo_scale_ratio, GM_Aur.AD_Leo_scale_ratio, Recx_11.AD_Leo_scale_ratio]
CTTS_RMSs = [TW_Hya.CTTS_RMS, GM_Aur.CTTS_RMS, Recx_11.CTTS_RMS]
WTTS_F_qs = [TW_Hya.F_q_wtts, GM_Aur.F_q_wtts, Recx_11.F_q_wtts]

inverse_CTTS_RMSs =[]
detectability = []
for i in range(len(scale_ratios)):
    detectability.append(WTTS_F_qs[i]/CTTS_RMSs[i])
    inverse_CTTS_RMSs.append(1/CTTS_RMSs[i])
    


num_of_5_sigma_expected_si4 = [TW_Hya.num_of_5_sigma_expected_si4, GM_Aur.num_of_5_sigma_expected_si4, Recx_11.num_of_5_sigma_expected_si4]
num_of_5_sigma_expected_si4_err_up =[TW_Hya.num_of_5_sigma_expected_si4_err_up, GM_Aur.num_of_5_sigma_expected_si4_err_up, Recx_11.num_of_5_sigma_expected_si4_err_up]
num_of_5_sigma_expected_si4_err_down = [TW_Hya.num_of_5_sigma_expected_si4_err_down, GM_Aur.num_of_5_sigma_expected_si4_err_down, Recx_11.num_of_5_sigma_expected_si4_err_down]

num_of_3_sigma_expected_si4 = [TW_Hya.num_of_3_sigma_expected_si4, GM_Aur.num_of_3_sigma_expected_si4, Recx_11.num_of_3_sigma_expected_si4]
num_of_3_sigma_expected_si4_err_up =[TW_Hya.num_of_3_sigma_expected_si4_err_up, GM_Aur.num_of_3_sigma_expected_si4_err_up, Recx_11.num_of_3_sigma_expected_si4_err_up]
num_of_3_sigma_expected_si4_err_down = [TW_Hya.num_of_3_sigma_expected_si4_err_down, GM_Aur.num_of_3_sigma_expected_si4_err_down, Recx_11.num_of_3_sigma_expected_si4_err_down]


plt.figure()
# plt.errorbar(Accretion_rates, num_of_5_sigma_expected_si4,
#              yerr = [num_of_5_sigma_expected_si4_err_down,num_of_5_sigma_expected_si4_err_up],
#              xerr = None, ecolor = 'red', capsize = 2, linewidth = 0, elinewidth = 1, marker = '.',
#              markersize = 6, color = 'green')
plt.errorbar(Accretion_rates, num_of_3_sigma_expected_si4,
             yerr = [num_of_3_sigma_expected_si4_err_down,num_of_3_sigma_expected_si4_err_up],
             xerr = None, ecolor = 'red', capsize = 2, linewidth = 0, elinewidth = 2, marker = '*',
             markersize = 9, color = 'green')

for i in range(len(Accretion_rates)):
    if i != 1:plt.text(Accretion_rates[i]+.2e-9,num_of_3_sigma_expected_si4[i], names[i], fontsize = 14, color = 'green')
    else: plt.text(Accretion_rates[i]-1.6e-9,num_of_3_sigma_expected_si4[i], names[i], fontsize = 14, color = 'green')

plt.ylabel('Flare Frequency  ($67 \/ ks^{-1}$)', fontsize = 18)
plt.xlabel ('Mass Accretion Rate ($M_{\u2609} \/ yr^{-1}$)', fontsize = 18)
plt.title('3$\sigma$ Flares in Si IV as a Function of Accretion Rate', fontsize = 18)




plt.figure()
plt.plot(Accretion_rates,detectability,marker = '.', linewidth = 0)


for i in range(len(Accretion_rates)):
    if i != 1:plt.text(Accretion_rates[i]+.2e-9,detectability[i], names[i], fontsize = 14, color = 'green')
    else: plt.text(Accretion_rates[i]-1.6e-9,detectability[i], names[i], fontsize = 14, color = 'green')

plt.ylabel('$\dfrac{F_{q,WTTS}}{RMS_{CTTS}}$', fontsize = 18)
plt.xlabel ('Mass Accretion Rate ($M_{\u2609} \/ yr^{-1}$)', fontsize = 18)
plt.title('detectability Estimator', fontsize = 18)



plt.figure()
plt.plot(Accretion_rates,inverse_CTTS_RMSs,marker = '.', linewidth = 0)


for i in range(len(Accretion_rates)):
    if i != 1:plt.text(Accretion_rates[i]+.2e-9,inverse_CTTS_RMSs[i], names[i], fontsize = 14, color = 'green')
    else: plt.text(Accretion_rates[i]-1.6e-9,inverse_CTTS_RMSs[i], names[i], fontsize = 14, color = 'green')

plt.ylabel('$(CTTS \/ \/ RMS)^{-1}$', fontsize = 18)
plt.xlabel ('Mass Accretion Rate ($M_{\u2609} \/ yr^{-1}$)', fontsize = 18)
plt.title('detectability Estimator', fontsize = 18)

plt.figure()
plt.errorbar(inverse_CTTS_RMSs, num_of_3_sigma_expected_si4,
             yerr = [num_of_3_sigma_expected_si4_err_down,num_of_3_sigma_expected_si4_err_up],
             xerr = None, ecolor = 'red', capsize = 2, linewidth = 0, elinewidth = 2, marker = '*',
             markersize = 9, color = 'green')

for i in range(len(Accretion_rates)):
    if i != 1:plt.text(inverse_CTTS_RMSs[i]+.2e-9,num_of_3_sigma_expected_si4[i], names[i], fontsize = 14, color = 'green')
    else: plt.text(inverse_CTTS_RMSs[i]-1.6e-9,num_of_3_sigma_expected_si4[i], names[i], fontsize = 14, color = 'green')

plt.ylabel('Flare Frequency  ($67 \/ ks^{-1}$)', fontsize = 18)
plt.xlabel ('$(CTTS \/ \/ RMS)^{-1}$', fontsize = 18)
plt.title('3$\sigma$ Flares in Si IV as a Function of Inverse RMS', fontsize = 18)

plt.figure()
plt.errorbar(CTTS_RMSs, num_of_3_sigma_expected_si4,
             yerr = [num_of_3_sigma_expected_si4_err_down,num_of_3_sigma_expected_si4_err_up],
             xerr = None, ecolor = 'red', capsize = 2, linewidth = 0, elinewidth = 2, marker = '*',
             markersize = 9, color = 'green')

for i in range(len(Accretion_rates)):
    if i != 1:plt.text(CTTS_RMSs[i]+.2e-9,num_of_3_sigma_expected_si4[i], names[i], fontsize = 14, color = 'green')
    else: plt.text(CTTS_RMSs[i]-1.6e-9,num_of_3_sigma_expected_si4[i], names[i], fontsize = 14, color = 'green')

plt.ylabel('Flare Frequency  ($67 \/ ks^{-1}$)', fontsize = 18)
plt.xlabel ('$(CTTS \/ \/ RMS)$', fontsize = 18)
plt.title('3$\sigma$ Flares in Si IV as a Function of RMS', fontsize = 18)







plt.close('all')
fig,ax = plt.subplots()
# make a plot
ax.errorbar(Accretion_rates, num_of_3_sigma_expected_si4,
             yerr = [num_of_3_sigma_expected_si4_err_down,num_of_3_sigma_expected_si4_err_up],
             xerr = None, ecolor = 'green', capsize = 7, linewidth = 0, elinewidth = 3, marker = '*',
             markersize = 13, color = 'green')

for i in range(len(Accretion_rates)):
    if i != 1:plt.text(Accretion_rates[i]+.2e-9,num_of_3_sigma_expected_si4[i], names[i], fontsize = 14, color = 'Black')
    else: plt.text(Accretion_rates[i]-1.6e-9,num_of_3_sigma_expected_si4[i], names[i], fontsize = 14, color = 'black')

plt.text(.6e-8,60,'$D \propto \dfrac{1}{CTTS_{RMS}}$', fontsize = 17, fontweight = 'bold')

ax.tick_params(axis='y', colors="green")

ax2=ax.twinx()
ax2.tick_params(axis='y', colors="blue")

# make a plot with different y-axis using second axis object
ax2.plot(Accretion_rates,detectability,marker = '*',markersize = 13, color = 'blue', linewidth = 0)

ax2.set_ylabel("Detectability",color="blue",fontsize=16)

# set x-axis label
ax.set_xlabel('Mass Accretion Rate ($M_{\u2609} \/ yr^{-1}$)',fontsize=16)
# set y-axis label
ax.set_ylabel("3$\sigma$ Si IV Flare Frequency  ($67 \/ ks^{-1}$)",color="green",fontsize=16)
ax.set_title('Flare Frequency and Detectability vs Accretion Rate',fontsize=16, fontweight = 'bold')









fig,ax = plt.subplots()
# make a plot
ax.errorbar(Accretion_rates, num_of_5_sigma_expected_si4,
             yerr = [num_of_5_sigma_expected_si4_err_down,num_of_5_sigma_expected_si4_err_up],
             xerr = None, ecolor = 'green', capsize = 7, linewidth = 0, elinewidth = 3, marker = '*',
             markersize = 13, color = 'green')

for i in range(len(Accretion_rates)):
    if i != 1:plt.text(Accretion_rates[i]+.2e-9,num_of_5_sigma_expected_si4[i], names[i], fontsize = 14, color = 'Black')
    else: plt.text(Accretion_rates[i]-1.6e-9,num_of_5_sigma_expected_si4[i], names[i], fontsize = 14, color = 'black')

plt.text(.6e-8,60,'$D \propto \dfrac{1}{CTTS_{RMS}}$', fontsize = 17, fontweight = 'bold')

ax.tick_params(axis='y', colors="green")

ax2=ax.twinx()
ax2.tick_params(axis='y', colors="blue")

# make a plot with different y-axis using second axis object
ax2.plot(Accretion_rates,detectability,marker = '*',markersize = 13, color = 'blue', linewidth = 0)

ax2.set_ylabel("Detectability",color="blue",fontsize=16)

# set x-axis label
ax.set_xlabel('Mass Accretion Rate ($M_{\u2609} \/ yr^{-1}$)',fontsize=16)
# set y-axis label
ax.set_ylabel("3$\sigma$ Si IV Flare Frequency  ($67 \/ ks^{-1}$)",color="green",fontsize=16)
ax.set_title('Flare Frequency and Detectability vs Accretion Rate',fontsize=16, fontweight = 'bold')





fig,ax = plt.subplots()
# make a plot
ax.errorbar(Accretion_rates, num_of_5_sigma_expected_si4,
             yerr = [num_of_5_sigma_expected_si4_err_down,num_of_5_sigma_expected_si4_err_up],
             xerr = None, ecolor = 'green', capsize = 7, linewidth = 0, elinewidth = 3, marker = '*',
             markersize = 13, color = 'green')

for i in range(len(Accretion_rates)):
    if i != 1:plt.text(Accretion_rates[i]+.2e-9,num_of_5_sigma_expected_si4[i], names[i], fontsize = 14, color = 'Black')
    else: plt.text(Accretion_rates[i]-1.6e-9,num_of_5_sigma_expected_si4[i], names[i], fontsize = 14, color = 'black')

plt.text(.6e-8,60,'$D \propto \dfrac{1}{CTTS_{RMS}}$', fontsize = 17, fontweight = 'bold')

ax.tick_params(axis='y', colors="green")

ax2=ax.twinx()
ax2.tick_params(axis='y', colors="red")

# make a plot with different y-axis using second axis object
ax2.plot(Accretion_rates,inverse_CTTS_RMSs,marker = '*',markersize = 13, color = 'red', linewidth = 0)

ax2.set_ylabel("Relative Detectability $(\sigma_{CTTS}^{-1})$",color="red",fontsize=16)

# set x-axis label
ax.set_xlabel('Mass Accretion Rate ($M_{\u2609} \/ yr^{-1}$)',fontsize=16)
# set y-axis label
ax.set_ylabel("Predicted # of 3$\sigma$ Si IV Flares $(67 \/ ks)^{-1}$",color="green",fontsize=16)
ax.set_title('Frequency and Detectability vs $\dot{M}$',fontsize=16)


