# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 10:31:44 2021

@author: pahi9557
"""


# -*- coding: utf-8 -*-
"""
NOTE! Any changes to line locations must be done in two places (total + subdivided x1d's')
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from tqdm import tqdm
import time as ostime
from scipy.optimize import curve_fit  


def func_linear(x, m, b):
    return (m*x+b)

def find_nearest(array,value):
    idx = np.nanargmin(np.abs(array - value))
    return idx

def onclick(event):
    global ix, iy, coords_x, coords_y, terminator
    ix, iy = event.xdata, event.ydata
    coords_x.append(ix)
    coords_y.append(iy)
    print('time:  '+ str(ix))
    if len(coords_x) == 2:
        fig.canvas.mpl_disconnect(fig)
        plt.close()
        coords_x = np.sort(coords_x)
        terminator = 1
    return

def subtract_quies(species_count_rate, quiescent_time_index):
    mean_quies = np.mean(species_count_rate[quiescent_time_index])
    subtracted_counts = species_count_rate - mean_quies
    return subtracted_counts

#function to find all tag files in directory
def get_file_names_with_strings(str_list):
    full_list = os.listdir(r'C:\Windows\System32\astroconda\AD_Leo_E140M\res_20s')
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]

    return final_list

#defining the time resolution of this file
time_resolution = 20

#finding and seperating all subdivided x1d files, and total x1d files
filename_list_all = get_file_names_with_strings(['x1d.fits']) #all of the files
filename_list_tot = get_file_names_with_strings(['_tot_']) #only the total x1d files
filename_list = list(set(filename_list_all) - set(filename_list_tot)) #all - tot = subdivided files
#now I want to sort them all according to their MJD, total files are already sorted

# This block of code was used to verify that the files I read in are in MJD
# sorted order

MJD_list = [] #defining empty MJD list, to be filled in 'for loop'
for i in filename_list:
    hdul = fits.open(r"C:\Windows\System32\astroconda\AD_Leo_E140M\res_20s\%s" % i)
    TEXPSTRT = hdul[0].header['TEXPSTRT'] #start time (MJD) of 1st exposure in file
    MJD_list.append(TEXPSTRT)
    hdul.close()
sort_key = np.argsort(MJD_list) # this is the order that filename_list should be in
temp = []   # reparing empty variable to be made into the new filename_list                             
for i in sort_key: temp.append(filename_list[i])
filename_list = temp                   
MJD_list = np.asarray(MJD_list); MJD_list = MJD_list[sort_key]
del temp; del filename_list_all; del sort_key; # getting rid of finsihed variables   
                
#preparing master lists
master_continuum_flux=np.asarray([]); master_continuum_counts=np.asarray([]); master_continuum_error=np.asarray([]);
master_si4_flux=np.asarray([]); master_si4_counts=np.asarray([]); master_si4_error=np.asarray([]);
master_c4_flux=np.asarray([]); master_c4_counts=np.asarray([]); master_c4_error=np.asarray([])
master_He2_flux=np.asarray([]); master_He2_counts=np.asarray([]); master_He2_error=np.asarray([]);



master_time=[]; time = 0
change_in_obs_flag = []       
                                  
file_count = 1
for i in filename_list: # this loop iterates through each sub-divided image, generates lightcurve
    print(' ');print(' ');
    print('Now working on %s of 26 E140M files' % file_count)
    change_in_obs_flag.append(time)
    hdul = fits.open(r"C:\Windows\System32\astroconda\AD_Leo_E140M\res_20s\%s" % i)
    HST_ID = i[0:9] #the first 9 characters of an HST filename is the ID
    TEXPTIME = hdul[0].header['TEXPTIME'] #total exposure time in seconds
    TEXPTIME = int(TEXPTIME) #rounding to the nearest second
    print('Total Exposure time =  %ss' % TEXPTIME)
    #locating the corresponding 'total' file for flux callibration
    filename_tot = '%s_raw_PH_tot_%ss_x1d.fits' % (HST_ID, TEXPTIME)
    try:
        hdul_tot  = fits.open(r"C:\Windows\System32\astroconda\AD_Leo_E140M\res_20s\%s" % filename_tot)
    except:
        print(''); print('Original attempt to open total file failed, conducting alternate attempt');
        ostime.sleep(2)
        filename_tot = get_file_names_with_strings([HST_ID + '_raw_PH_tot_'])
        filename_tot = filename_tot[0]
        hdul_tot  = fits.open(r"C:\Windows\System32\astroconda\AD_Leo_E140M\res_20s\%s" % filename_tot)
    #first to extract the data from the total file
    
    tot_x1d_wav_array = []
    tot_x1d_flux_array = []
    tot_x1d_count_array = []
    
    #extracting average count rates and fluxes from total x1d
    #this same extraction method will be looped over afterwards for sub x1d exposures
    data = hdul_tot['sci',1].data
    tot_x1d_wav_array.append(data['wavelength'])
    tot_x1d_flux_array.append(data['flux'])
    tot_x1d_count_array.append(data['net'])
    #this previous step may seem unnecessary, but this is what is required later down
    #the line so I just kept the same style here
    wave = tot_x1d_wav_array[0]; wave = wave.flatten()
    flux = tot_x1d_flux_array[0]; flux = flux.flatten()
    counts = tot_x1d_count_array[0]; counts = counts.flatten()
    
    sort_key = np.argsort(wave)
    wave = wave[sort_key]; flux = flux[sort_key]; counts = counts[sort_key]
    
    ################ indentifying line locations############
    c3_line_location = np.where((wave > 1174) & (wave < 1176.5))
    n5_line_location = np.where((wave > 1238) & (wave < 1243.3))
    c4_line_location = np.where((wave > 1547.6) & (wave < 1551.5))
    He2_line_location = np.where((wave > 1638) & (wave < 1642))
    
    c2_line_location = np.where((wave > 1334) & (wave < 1336.5))
    si2_line_location = np.where((wave > 1264.34) & (wave < 1265.5))
    si3_line_location = np.where((wave > 1205.9) & (wave < 1207.3))
    o1_line_location = np.where((wave > 1304.6) & (wave < 1306.5))
    
    continuum_line_location = np.where((wave > 1339) & (wave < 1351))
    
    bandpass1_1 =  np.asarray(np.where((wave > 1170 ) & (wave < 1210))) #stop at Ly Alpha
    bandpass1_2 =  np.asarray(np.where((wave > 1220 ) & (wave < 1300))) #stop at O I airglow
    bandpass1_3 =  np.asarray(np.where((wave > 1310 ) & (wave < 1410))) #stop at bandpass2
    #np.where makes lame tuples so I have to use this funny method to append them all together
    bandpass1_line_location = []
    for j in range(np.size(bandpass1_1)):
        bandpass1_line_location.append(bandpass1_1[0][j])
    for j in range(np.size(bandpass1_2)):
        bandpass1_line_location.append(bandpass1_2[0][j])
    for j in range(np.size(bandpass1_3)):
        bandpass1_line_location.append(bandpass1_3[0][j])
    bandpass2_line_location = np.where((wave > 1410 ) & (wave < 1680))

    #now for the doublets
    si4_line_location_1 = np.where((wave > 1393) & (wave < 1395))
    si4_line_location_2 = np.where((wave > 1402) & (wave < 1404)) 
    si4_line_location = []
    for j in range(np.size(si4_line_location_1)):
        si4_line_location.append(si4_line_location_1[0][j])
    for j in range(np.size(si4_line_location_2)):
        si4_line_location.append(si4_line_location_2[0][j])
    del si4_line_location_1; del si4_line_location_2;
    ##############################################################
    
    
    #Now to compute the average fluxes and count rates from total x1d file#
    
    #avg flux
    avg_c2_line_flux = np.trapz(flux[c2_line_location],wave[c2_line_location])
    avg_c3_line_flux = np.trapz(flux[c3_line_location],wave[c3_line_location])
    avg_n5_line_flux = np.trapz(flux[n5_line_location],wave[n5_line_location])
    avg_si4_line_flux = np.trapz(flux[si4_line_location],wave[si4_line_location])
    avg_c4_line_flux = np.trapz(flux[c4_line_location],wave[c4_line_location])
    avg_He2_line_flux = np.trapz(flux[He2_line_location],wave[He2_line_location])
    avg_continuum_line_flux = np.trapz(flux[continuum_line_location],wave[continuum_line_location])
    avg_bandpass1_line_flux = np.trapz(flux[bandpass1_line_location],wave[bandpass1_line_location])
    avg_bandpass2_line_flux = np.trapz(flux[bandpass2_line_location],wave[bandpass2_line_location])   
    #avg counts
    avg_c2_line_counts = np.sum(counts[c2_line_location])
    avg_c3_line_counts = np.sum(counts[c3_line_location])
    avg_n5_line_counts = np.sum(counts[n5_line_location])
    avg_si4_line_counts = np.sum(counts[si4_line_location])
    avg_c4_line_counts = np.sum(counts[c4_line_location])
    avg_He2_line_counts = np.sum(counts[He2_line_location])
    avg_continuum_line_counts = np.sum(counts[continuum_line_location])
    avg_bandpass1_line_counts = np.sum(counts[bandpass1_line_location])
    avg_bandpass2_line_counts = np.sum(counts[bandpass2_line_location])   
    
    #Now that I have computed the average flux and average counts at each line for the
    #total exposure, I want to calculate the values at each time step so that I can
    #generate light curves
    
    number_of_images = hdul[0].header['NEXTEND'] # number of frames
    x1d_wav_array = []
    x1d_count_array = []
    x1d_error_array = []
    
    #extracting the useful data from the fits file
    for i in range(number_of_images):
        data = hdul['sci',i+1].data
        x1d_wav_array.append(data['wavelength'])
        x1d_count_array.append(data['net'])
    
    #preparing empty variables for each line
    c3_flux=[]; c3_flux_error=[]; c3_counts=[]; c3_error=[];
    n5_flux=[]; n5_flux_error=[]; n5_counts=[]; n5_error=[]
    si4_flux=[]; si4_flux_error=[]; si4_counts=[]; si4_error=[]
    c4_flux=[]; c4_flux_error=[]; c4_counts=[]; c4_error=[]
    He2_flux=[]; He2_flux_error=[]; He2_counts=[]; He2_error=[]
    #adding extra lines to compare with Hawley et al 2003
    c2_flux=[]; c2_flux_error=[]; c2_counts=[]; c2_error=[];
    si2_flux=[]; si2_flux_error=[]; si2_counts=[]; si2_error=[];
    si3_flux=[]; si3_flux_error=[]; si3_counts=[]; si3_error=[];
    o1_flux=[]; o1_flux_error=[]; o1_counts=[]; o1_error=[];
    #larger non-line regions
    continuum_flux=[]; continuum_flux_error=[]; continuum_counts=[]; continuum_error=[];
    bandpass1_flux=[]; bandpass1_flux_error=[]; bandpass1_counts=[]; bandpass1_error=[];
    bandpass2_flux=[]; bandpass2_flux_error=[]; bandpass2_counts=[]; bandpass2_error=[];
    
    

    
    #for each sub exposure: computing flux, err, and count rate
    for i in tqdm(range(len(x1d_wav_array))):
        
        wave = x1d_wav_array[i]; wave = wave.flatten()
        counts = x1d_count_array[i]; counts = counts.flatten()
    
        sort_key = np.argsort(wave)
        wave = wave[sort_key]; counts = counts[sort_key]
    
    
        ################ indentifying line locations############
        c3_line_location = np.where((wave > 1174) & (wave < 1176.5))
        n5_line_location = np.where((wave > 1238) & (wave < 1243.3))
        c4_line_location = np.where((wave > 1547.6) & (wave < 1551.5))
        He2_line_location = np.where((wave > 1639.6) & (wave < 1641.5))
        
        c2_line_location = np.where((wave > 1334) & (wave < 1336.5))
        si2_line_location = np.where((wave > 1264.34) & (wave < 1265.5))
        si3_line_location = np.where((wave > 1205.9) & (wave < 1207.3))
        o1_line_location = np.where((wave > 1304.6) & (wave < 1306.5))
        
        continuum_line_location = np.where((wave > 1339) & (wave < 1351))
        
        bandpass1_1 =  np.asarray(np.where((wave > 1170 ) & (wave < 1210))) #stop at Ly Alpha
        bandpass1_2 =  np.asarray(np.where((wave > 1220 ) & (wave < 1300))) #stop at O I airglow
        bandpass1_3 =  np.asarray(np.where((wave > 1310 ) & (wave < 1410))) #stop at bandpass2
        #np.where makes lame tuples so I have to use this funny method to append them all together
        bandpass1_line_location = []
        for j in range(np.size(bandpass1_1)):
            bandpass1_line_location.append(bandpass1_1[0][j])
        for j in range(np.size(bandpass1_2)):
            bandpass1_line_location.append(bandpass1_2[0][j])
        for j in range(np.size(bandpass1_3)):
            bandpass1_line_location.append(bandpass1_3[0][j])
        bandpass2_line_location = np.where((wave > 1410 ) & (wave < 1680))
    
        #now for the doublets
        si4_line_location_1 = np.where((wave > 1393) & (wave < 1395))
        si4_line_location_2 = np.where((wave > 1402) & (wave < 1404)) 
        si4_line_location = []
        for j in range(np.size(si4_line_location_1)):
            si4_line_location.append(si4_line_location_1[0][j])
        for j in range(np.size(si4_line_location_2)):
            si4_line_location.append(si4_line_location_2[0][j])
        del si4_line_location_1; del si4_line_location_2;
        ################################################################
        
        
        
        # computing count rates
        c2temp_counts = np.sum(counts[c2_line_location])
        c3temp_counts = np.sum(counts[c3_line_location])
        c4temp_counts = np.sum(counts[c4_line_location])
        
        continuum_temp_counts = np.sum(counts[continuum_line_location])
        si4_temp_counts = np.sum(counts[si4_line_location])
        He2_temp_counts = np.sum(counts[He2_line_location])
        
        #appending count rates
        c2_counts.append(c2temp_counts); c2_error.append(np.sqrt(c2temp_counts*time_resolution)/time_resolution)
        c3_counts.append(c3temp_counts); c3_error.append(np.sqrt(c3temp_counts*time_resolution)/time_resolution)
        c4_counts.append(c4temp_counts); c4_error.append(np.sqrt(c4temp_counts*time_resolution)/time_resolution)
        
        continuum_counts.append(continuum_temp_counts); continuum_error.append(np.sqrt(continuum_temp_counts*time_resolution)/time_resolution)
        si4_counts.append(si4_temp_counts); si4_error.append(np.sqrt(si4_temp_counts*time_resolution)/time_resolution)
        He2_counts.append(He2_temp_counts); He2_error.append(np.sqrt(He2_temp_counts*time_resolution)/time_resolution)
        
        
        
        master_time.append(time)
        time = time + time_resolution
            
     #converting count lists to arrays
    c2_counts=np.asarray(c2_counts); c2_error=np.asarray(c2_error);
    c3_counts=np.asarray(c3_counts); c3_error=np.asarray(c3_error);
    c4_counts=np.asarray(c4_counts); c4_error=np.asarray(c4_error);
    
    continuum_counts=np.asarray(continuum_counts); continuum_error=np.asarray(continuum_error);
    si4_counts=np.asarray(si4_counts); si4_error=np.asarray(si4_error);
    He2_counts=np.asarray(He2_counts); He2_error=np.asarray(He2_error);
    
    # Now that we have the average flux and avg count rate, as well as the count
    #rate in each sub exposure: we can scale the sub exposure rates by the average 
    #rates to effectivley flux-callibrate our curves
    
    fluxed_c2 = c2_counts*(avg_c2_line_flux/avg_c2_line_counts)
    fluxed_c2_err = c2_error*(avg_c2_line_flux/avg_c2_line_counts)
    fluxed_c3 = c3_counts*(avg_c3_line_flux/avg_c3_line_counts)
    fluxed_c3_err = c3_error*(avg_c3_line_flux/avg_c3_line_counts)
    fluxed_c4 = c4_counts*(avg_c4_line_flux/avg_c4_line_counts)
    fluxed_c4_err = c4_error*(avg_c4_line_flux/avg_c4_line_counts)


    fluxed_continuum = continuum_counts*(avg_continuum_line_flux/avg_continuum_line_counts)
    fluxed_continuum_err = continuum_error*(avg_continuum_line_flux/avg_continuum_line_counts)
    fluxed_si4 = si4_counts*(avg_si4_line_flux/avg_si4_line_counts)
    fluxed_si4_err = si4_error*(avg_si4_line_flux/avg_si4_line_counts)
    fluxed_He2 = He2_counts*(avg_He2_line_flux/avg_He2_line_counts)
    fluxed_He2_err = He2_error*(avg_He2_line_flux/avg_He2_line_counts)


    master_continuum_counts= np.append(master_continuum_counts,continuum_counts);
    master_continuum_flux=   np.append(master_continuum_flux,fluxed_continuum); 
    master_continuum_error=  np.append(master_continuum_error,fluxed_continuum_err); 
    master_si4_counts= np.append(master_si4_counts,si4_counts);
    master_si4_flux=   np.append(master_si4_flux,fluxed_si4); 
    master_si4_error=  np.append(master_si4_error,fluxed_si4_err); 
    
    master_c4_counts= np.append(master_c4_counts,c4_counts);
    master_c4_flux=   np.append(master_c4_flux,fluxed_c4); 
    master_c4_error=  np.append(master_c4_error,fluxed_c4_err); 
    
    master_He2_counts= np.append(master_He2_counts,He2_counts);
    master_He2_flux=   np.append(master_He2_flux,fluxed_He2); 
    master_He2_error=  np.append(master_He2_error,fluxed_He2_err); 




    file_count += 1  
    
#deleting variables used in loops
del He2_counts;del He2_error;del He2_flux;del He2_flux_error;del He2_line_location;
del avg_He2_line_counts;del avg_He2_line_flux;del avg_bandpass1_line_counts;
del avg_bandpass1_line_flux;del avg_bandpass2_line_counts;del avg_bandpass2_line_flux;
del avg_c2_line_flux;del avg_c3_line_counts;del avg_c3_line_flux;del avg_c4_line_counts;
del avg_c4_line_flux;del avg_continuum_line_counts;del avg_continuum_line_flux;del avg_n5_line_counts;
del avg_n5_line_flux;del avg_si4_line_counts;del avg_si4_line_flux;del bandpass1_1;del bandpass1_2;
del bandpass1_3;del bandpass1_counts;del bandpass1_error;del bandpass1_flux;del bandpass1_flux_error;
del bandpass1_line_location;del bandpass2_counts;del bandpass2_error;del bandpass2_flux;del bandpass2_flux_error;
del bandpass2_line_location;del c2_counts;del c2_error;del c2_flux;del c2_flux_error;del c2_line_location;
del c2temp_counts;del c3_counts;del c3_error;del c3_flux;del i;del j;del avg_c2_line_counts;del c3_flux_error;
del c3_line_location;del c3temp_counts;
    

#scaling
ADLEO_TWA7_si4_flux = master_si4_flux*0.15*0.321   #*0.321 is flux distance correction
ADLEO_TWA7_si4_error = master_si4_error*0.15*0.321  

ADLEO_TWA7_c4_flux = master_c4_flux*0.075*0.321      
ADLEO_TWA7_c4_error = master_c4_error*0.075*0.321   

ADLEO_TWA7_He2_flux = master_He2_flux*0.13*0.321      
ADLEO_TWA7_He2_error = master_He2_error*0.13*0.321   

master_time = np.asarray(master_time); change_in_obs_flag = np.asarray(change_in_obs_flag);
master_time = master_time/1000; change_in_obs_flag = change_in_obs_flag/1000



##############################################
#whole light curve plot with all of the species, no CTTS lines y axis 0:1
fig,subplot = plt.subplots(3, 1)
filename_counter = 0

# for i in change_in_obs_flag:
    # flag = i in subplot0_HST_ID_fixes_list
    # if i < master_time[1117] and flag == False:
    #     subplot[0].axvline(x=i, color = 'green', alpha = 0.25, linewidth = 2)
    #     subplot[0].text(x=i,y=0.65,s=filename_list[filename_counter][0:9], color = 'green', alpha = 0.75)
    #     subplot[0].text(x=i,y=0.58,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75)
    # if i < master_time[1117]: filename_counter+=1
subplot[0].plot(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117]-ADLEO_TWA7_si4_error[0:1117],ADLEO_TWA7_si4_flux[0:1117]+ ADLEO_TWA7_si4_error[0:1117],color = 'green', alpha = 0.1)

subplot[0].plot(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]-ADLEO_TWA7_c4_error[0:1117],ADLEO_TWA7_c4_flux[0:1117]+ ADLEO_TWA7_c4_error[0:1117],color = 'blue', alpha = 0.1)

subplot[0].plot(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117]-ADLEO_TWA7_He2_error[0:1117],ADLEO_TWA7_He2_flux[0:1117]+ ADLEO_TWA7_He2_error[0:1117],color = 'red', alpha = 0.1)

subplot[0].legend(fontsize = 13, loc = 'upper left')
# subplot[0].set_ylim(0,1)



subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234]-ADLEO_TWA7_si4_error[1117:2234],ADLEO_TWA7_si4_flux[1117:2234]+ ADLEO_TWA7_si4_error[1117:2234],color = 'green', alpha = 0.1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]-ADLEO_TWA7_c4_error[1117:2234],ADLEO_TWA7_c4_flux[1117:2234]+ ADLEO_TWA7_c4_error[1117:2234],color = 'blue', alpha = 0.1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234]-ADLEO_TWA7_He2_error[1117:2234],ADLEO_TWA7_He2_flux[1117:2234]+ ADLEO_TWA7_He2_error[1117:2234],color = 'red', alpha = 0.1)



subplot[2].plot(master_time[2234:], ADLEO_TWA7_si4_flux[2234:], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_si4_flux[2234:]-ADLEO_TWA7_si4_error[2234:],ADLEO_TWA7_si4_flux[2234:]+ ADLEO_TWA7_si4_error[2234:],color = 'green', alpha = 0.1)

subplot[2].plot(master_time[2234:], ADLEO_TWA7_c4_flux[2234:], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]-ADLEO_TWA7_c4_error[2234:],ADLEO_TWA7_c4_flux[2234:]+ ADLEO_TWA7_c4_error[2234:],color = 'blue', alpha = 0.1)

subplot[2].plot(master_time[2234:], ADLEO_TWA7_He2_flux[2234:], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_He2_flux[2234:]-ADLEO_TWA7_He2_error[2234:],ADLEO_TWA7_He2_flux[2234:]+ ADLEO_TWA7_He2_error[2234:],color = 'red', alpha = 0.1)


# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Flux  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to quiescent level of TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 

##############################################
#whole light curve plot with all of the species, no CTTS lines, subtract the quiescence flux
fig,subplot = plt.subplots(3, 1)
filename_counter = 0

# for i in change_in_obs_flag:
    # flag = i in subplot0_HST_ID_fixes_list
    # if i < master_time[1117] and flag == False:
    #     subplot[0].axvline(x=i, color = 'green', alpha = 0.25, linewidth = 2)
    #     subplot[0].text(x=i,y=0.65,s=filename_list[filename_counter][0:9], color = 'green', alpha = 0.75)
    #     subplot[0].text(x=i,y=0.58,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75)
    # if i < master_time[1117]: filename_counter+=1
subplot[0].plot(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117]-1.5e-14, color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117]-ADLEO_TWA7_si4_error[0:1117]-1.5e-14,ADLEO_TWA7_si4_flux[0:1117]+ ADLEO_TWA7_si4_error[0:1117]-1.5e-14,color = 'green', alpha = 0.1)

subplot[0].plot(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]-3e-14, color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]-ADLEO_TWA7_c4_error[0:1117]-3e-14,ADLEO_TWA7_c4_flux[0:1117]+ ADLEO_TWA7_c4_error[0:1117]-3e-14,color = 'blue', alpha = 0.1)

subplot[0].plot(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117]-3.2e-14, color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117]-ADLEO_TWA7_He2_error[0:1117]-3.2e-14,ADLEO_TWA7_He2_flux[0:1117]+ ADLEO_TWA7_He2_error[0:1117]-3.2e-14,color = 'red', alpha = 0.1)

subplot[0].axhline(y=0, color = 'black', alpha = 0.25, linewidth = 2, linestyle = 'dashed')

subplot[0].legend(fontsize = 13, loc = 'upper left')
# subplot[0].set_ylim(0,1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234]-1.5e-14, color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234]-ADLEO_TWA7_si4_error[1117:2234]-1.5e-14,ADLEO_TWA7_si4_flux[1117:2234]+ ADLEO_TWA7_si4_error[1117:2234]-1.5e-14,color = 'green', alpha = 0.1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]-3e-14, color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]-ADLEO_TWA7_c4_error[1117:2234]-3e-14,ADLEO_TWA7_c4_flux[1117:2234]+ ADLEO_TWA7_c4_error[1117:2234]-3e-14,color = 'blue', alpha = 0.1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234]-3.2e-14, color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234]-ADLEO_TWA7_He2_error[1117:2234]-3.2e-14,ADLEO_TWA7_He2_flux[1117:2234]+ ADLEO_TWA7_He2_error[1117:2234]-3.2e-14,color = 'red', alpha = 0.1)

subplot[1].axhline(y=0, color = 'black', alpha = 0.25, linewidth = 2, linestyle = 'dashed')

subplot[2].plot(master_time[2234:], ADLEO_TWA7_si4_flux[2234:]-1.5e-14, color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_si4_flux[2234:]-ADLEO_TWA7_si4_error[2234:]-1.5e-14,ADLEO_TWA7_si4_flux[2234:]+ ADLEO_TWA7_si4_error[2234:]-1.5e-14,color = 'green', alpha = 0.1)

subplot[2].plot(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]-3e-14, color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]-ADLEO_TWA7_c4_error[2234:]-3e-14,ADLEO_TWA7_c4_flux[2234:]+ ADLEO_TWA7_c4_error[2234:]-3e-14,color = 'blue', alpha = 0.1)

subplot[2].plot(master_time[2234:], ADLEO_TWA7_He2_flux[2234:]-3.2e-14, color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_He2_flux[2234:]-ADLEO_TWA7_He2_error[2234:]-3.2e-14,ADLEO_TWA7_He2_flux[2234:]+ ADLEO_TWA7_He2_error[2234:]-3.2e-14,color = 'red', alpha = 0.1)


subplot[2].axhline(y=0, color = 'black', alpha = 0.25, linewidth = 2, linestyle = 'dashed')


# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'$Flux - F_q$  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to quiescent level of TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 



##############################################
#whole light curve plot with all of the species, y axis 0:1
fig,subplot = plt.subplots(3, 1)
filename_counter = 0

# for i in change_in_obs_flag:
    # flag = i in subplot0_HST_ID_fixes_list
    # if i < master_time[1117] and flag == False:
    #     subplot[0].axvline(x=i, color = 'green', alpha = 0.25, linewidth = 2)
    #     subplot[0].text(x=i,y=0.65,s=filename_list[filename_counter][0:9], color = 'green', alpha = 0.75)
    #     subplot[0].text(x=i,y=0.58,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75)
    # if i < master_time[1117]: filename_counter+=1
subplot[0].plot(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117]-ADLEO_TWA7_si4_error[0:1117],ADLEO_TWA7_si4_flux[0:1117]+ ADLEO_TWA7_si4_error[0:1117],color = 'green', alpha = 0.1)

subplot[0].plot(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]-ADLEO_TWA7_c4_error[0:1117],ADLEO_TWA7_c4_flux[0:1117]+ ADLEO_TWA7_c4_error[0:1117],color = 'blue', alpha = 0.1)

subplot[0].plot(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117]-ADLEO_TWA7_He2_error[0:1117],ADLEO_TWA7_He2_flux[0:1117]+ ADLEO_TWA7_He2_error[0:1117],color = 'red', alpha = 0.1)

subplot[0].legend(fontsize = 13, loc = 'upper left')
# subplot[0].set_ylim(0,1)



subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234]-ADLEO_TWA7_si4_error[1117:2234],ADLEO_TWA7_si4_flux[1117:2234]+ ADLEO_TWA7_si4_error[1117:2234],color = 'green', alpha = 0.1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]-ADLEO_TWA7_c4_error[1117:2234],ADLEO_TWA7_c4_flux[1117:2234]+ ADLEO_TWA7_c4_error[1117:2234],color = 'blue', alpha = 0.1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234]-ADLEO_TWA7_He2_error[1117:2234],ADLEO_TWA7_He2_flux[1117:2234]+ ADLEO_TWA7_He2_error[1117:2234],color = 'red', alpha = 0.1)



subplot[2].plot(master_time[2234:], ADLEO_TWA7_si4_flux[2234:], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_si4_flux[2234:]-ADLEO_TWA7_si4_error[2234:],ADLEO_TWA7_si4_flux[2234:]+ ADLEO_TWA7_si4_error[2234:],color = 'green', alpha = 0.1)

subplot[2].plot(master_time[2234:], ADLEO_TWA7_c4_flux[2234:], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]-ADLEO_TWA7_c4_error[2234:],ADLEO_TWA7_c4_flux[2234:]+ ADLEO_TWA7_c4_error[2234:],color = 'blue', alpha = 0.1)

subplot[2].plot(master_time[2234:], ADLEO_TWA7_He2_flux[2234:], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_He2_flux[2234:]-ADLEO_TWA7_He2_error[2234:],ADLEO_TWA7_He2_flux[2234:]+ ADLEO_TWA7_He2_error[2234:],color = 'red', alpha = 0.1)

#horizantal lines for CTTS quiescent levels

#lows
subplot[0].axhline(y=0.7e-12, color = 'green', alpha = 0.25, linewidth = 2, linestyle = 'dashed')
subplot[0].axhline(y=1.5e-12, color = 'blue', alpha = 0.25, linewidth = 2, linestyle = 'dashed')
subplot[0].axhline(y=0.65e-12, color = 'red', alpha = 0.25, linewidth = 2, linestyle = 'dashed')

subplot[1].axhline(y=0.7e-12, color = 'green', alpha = 0.25, linewidth = 2, linestyle = 'dashed')
subplot[1].axhline(y=1.5e-12, color = 'blue', alpha = 0.25, linewidth = 2, linestyle = 'dashed')
subplot[1].axhline(y=0.65e-12, color = 'red', alpha = 0.25, linewidth = 2, linestyle = 'dashed')

subplot[2].axhline(y=0.7e-12, color = 'green', alpha = 0.25, linewidth = 2, linestyle = 'dashed')
subplot[2].axhline(y=1.5e-12, color = 'blue', alpha = 0.25, linewidth = 2, linestyle = 'dashed')
subplot[2].axhline(y=0.65e-12, color = 'red', alpha = 0.25, linewidth = 2, linestyle = 'dashed')

#highs
subplot[0].axhline(y=1.3e-12, color = 'green', alpha = 0.25, linewidth = 2)
subplot[0].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)
subplot[0].axhline(y=1.7e-12, color = 'red', alpha = 0.25, linewidth = 2)

subplot[1].axhline(y=1.3e-12, color = 'green', alpha = 0.25, linewidth = 2)
subplot[1].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)
subplot[1].axhline(y=1.7e-12, color = 'red', alpha = 0.25, linewidth = 2)

subplot[2].axhline(y=1.3e-12, color = 'green', alpha = 0.25, linewidth = 2)
subplot[2].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)
subplot[2].axhline(y=1.7e-12, color = 'red', alpha = 0.25, linewidth = 2)





# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Flux  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 




###plotting just cIV


##############################################
#whole light curve plot with all of the species, y axis 0:1
fig,subplot = plt.subplots(3, 1)
filename_counter = 0


subplot[0].plot(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]-ADLEO_TWA7_c4_error[0:1117],ADLEO_TWA7_c4_flux[0:1117]+ ADLEO_TWA7_c4_error[0:1117],color = 'blue', alpha = 0.1)

subplot[0].legend(fontsize = 13, loc = 'upper left')
# subplot[0].set_ylim(0,1)



subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]-ADLEO_TWA7_c4_error[1117:2234],ADLEO_TWA7_c4_flux[1117:2234]+ ADLEO_TWA7_c4_error[1117:2234],color = 'blue', alpha = 0.1)


subplot[2].plot(master_time[2234:], ADLEO_TWA7_c4_flux[2234:], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]-ADLEO_TWA7_c4_error[2234:],ADLEO_TWA7_c4_flux[2234:]+ ADLEO_TWA7_c4_error[2234:],color = 'blue', alpha = 0.1)

#horizantal lines for CTTS quiescent levels

#lows

subplot[0].axhline(y=1.5e-12, color = 'blue', alpha = 0.25, linewidth = 2, linestyle = 'dashed')

subplot[1].axhline(y=1.5e-12, color = 'blue', alpha = 0.25, linewidth = 2, linestyle = 'dashed')

subplot[2].axhline(y=1.5e-12, color = 'blue', alpha = 0.25, linewidth = 2, linestyle = 'dashed')

#highs
subplot[0].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)


subplot[1].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)

subplot[2].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)


# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Flux  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 




#flare locations from AD_Leo_FFD.py
peaks=np.array([  55, 111,  134,  246,  264,  289,  369,  564,  694,  718,  783,  808,
        841,  902, 1066, 1128, 1178, 1226, 1249, 1293, 1443, 1496, 1547,
       1566, 1578, 1608, 1630, 1716, 1774, 1816, 1892, 1909, 1989, 2006,
       2023, 2046, 2094, 2268, 2349, 2480, 2522, 2549, 2591, 2617, 2730,
       2769, 2806, 2823, 2876, 2908, 3019, 3049, 3110, 3145, 3164]) #added 111 and
      #potentially some others, so cross compare with AD Leo before altering


peaks_c4= peaks
peaks_si4 = np.array([  55,   99,  167,  246,  264,  289,  303,  355,  369,  398,  430,
        516,  565,  595,  695,  718,  783,  808,  839,  853,  904,  961,
       1064, 1128, 1181, 1226, 1248, 1292, 1337, 1443, 1496, 1524, 1546,
       1564, 1578, 1592, 1608, 1629, 1716, 1774, 1815, 1850, 1863, 1890,
       1906, 1934, 1950, 1989, 2003, 2023, 2046, 2064, 2094, 2194, 2234,
       2268, 2345, 2442, 2478, 2495, 2523, 2549, 2591, 2617, 2670, 2766,
       2823, 2908, 3019, 3033, 3049, 3110, 3148, 3163, 3283, 3318])
peaks_He2 = np.array([ 157,  262,  282,  782,  821,  906, 1226, 1248, 1500, 1578, 1610,
       1775, 1891, 1913, 1989, 2008, 2049, 2522, 2677, 2769, 3021])

#min/max limiting scenario
popt_si4_low = (1*10**(-12))*np.array([0, 0, 0.0139226519]) # computed in gaussian stats file
popt_si4_high = (1*10**(-12))*np.array([0,0 , 0.0444758789]) # computed in gaussian stats file

popt_c4_low = (1*10**(-12))*np.array([6.74862555, 1.535871794 , 0.048451233]) # computed in gaussian stats file
popt_c4_high = (1*10**(-12))*np.array([7.34673491, 4.0389487 , 0.06856082]) # computed in gaussian stats file

popt_He2_low = (1*10**(-12))*np.array([0,0, 0.0284057926]) # computed in gaussian stats file
popt_He2_high = (1*10**(-12))*np.array([0,0, 0.05246809048]) # computed in gaussian stats file

med =np.median(ADLEO_TWA7_c4_flux)

print('\n');print('\n');print('\n')
print('now computing estimations for the number of flares we expect to see above' 
      'the CTTS RMS when AD Leo is scaled to the level of the WTTS and and added'
       'to the CTTS signal. \n !!!WARNING!!! \n these are old results, new ones have'
       'been computed')

#limiting case for si4!!
ADLEO_TWA7_si4_quies = 1.7*10**-14
#after multiplying by 10^14 to bring the the flux up to normal magnitudes (0-3)
#the quiesence is noe 1.7, so we just divide all of our values by this to
#  "normalize" w.r.t. 1


observed_data = (ADLEO_TWA7_si4_flux-1.5e-14)*10**14 + 1 # add one b/c we're normalizing eveyrthing to 1
#and it's currently normed to 0
observed_data_err = (ADLEO_TWA7_si4_error)*10**14
popt_si4_low = (popt_si4_low)*10**14
popt_si4_high = (popt_si4_high)*10**14

# observed_data = (ADLEO_TWA7_si4_flux/1.7)*10**14 
# observed_data_err = (ADLEO_TWA7_si4_error/1.7)*10**14
# popt_si4_low = (popt_si4_low/1.7)*10**14
# popt_si4_high = (popt_si4_high/1.7)*10**14



fig,subplot = plt.subplots(2, 1)

subplot[0].errorbar(master_time, observed_data, yerr = np.sqrt(observed_data_err**2+popt_si4_low[2]**2),
             ecolor = 'grey', color = 'black', elinewidth = 0.5, linewidth = 0,
             capsize = 3, marker = '.', markersize = 3, label = r'Si IV (1394 + 1403)$\AA$')

subplot[0].plot(master_time,
                np.ones(len(master_time))*3*popt_si4_low[2]+1,
                color = 'orange', linewidth = 1, label = r'3$\sigma$', linestyle = 'dashed')
subplot[0].plot(master_time,
                np.ones(len(master_time))*5*popt_si4_low[2]+1,
                color = 'red', linewidth = 1, label = r'5$\sigma$')
subplot[0].plot(master_time[peaks_si4], observed_data[peaks_si4], "xr")

subplot[0].text(60,50,'low rms case', fontsize = 17, fontstyle = 'italic')
subplot[0].legend(fontsize = 13, loc = 'upper left')


subplot[1].errorbar(master_time, observed_data, yerr =  np.sqrt(observed_data_err**2+popt_si4_high[2]**2),
             ecolor = 'grey', color = 'black', elinewidth = 0.5, linewidth = 0,
             capsize = 3, marker = '.', markersize = 3, label = r'C IV (1548.2 + 1550.7)$\AA$')

subplot[1].plot(master_time,
                np.ones(len(master_time))*3*popt_si4_high[2]+1,
                color = 'orange', linewidth = 1, label = r'3$\sigma$ high', linestyle = 'dashed')
subplot[1].plot(master_time,
                np.ones(len(master_time))*5*popt_si4_high[2]+1,
                color = 'red', linewidth = 1, label = r'5$\sigma$ high')
subplot[1].plot(master_time[peaks_si4], observed_data[peaks_si4], "xr")


subplot[1].text(60,50,'high rms case', fontsize = 17, fontstyle = 'italic')

fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Scaled Flux', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'TW Hya Si IV Flare Propensity ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)



three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_si4)):
    temp = observed_data[peaks_si4[i]]
    if temp > 3*popt_si4_high[2]+1:
        three_sigma_list.append(i)
    if temp > 5*popt_si4_high[2]+1:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_si4[i]] + np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_high[2]**2)
    if temp > 3*popt_si4_high[2]+1:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_si4_high[2]+1:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_si4[i]] - np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_high[2]**2)
    if temp > 3*popt_si4_high[2]+1:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_si4_high[2]+1:
        five_sigma_list_err_down.append(i)

print("for Si IV high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))

three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_si4)):
    temp = observed_data[peaks_si4[i]]
    if temp > 3*popt_si4_low[2]+1:
        three_sigma_list.append(i)
    if temp > 5*popt_si4_low[2]+1:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_si4[i]] + np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_low[2]**2)
    if temp > 3*popt_si4_low[2]+1:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_si4_low[2]+1:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_si4[i]] - np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_low[2]**2)
    if temp > 3*popt_si4_low[2]+1:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_si4_low[2]+1:
        five_sigma_list_err_down.append(i)

print("for Si IV low:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))









#V2 of limiting case for c4

ADLEO_TWA7_c4_quies = 3.6*10**-14
#after multiplying by 10^14 to bring the the flux up to normal magnitudes (0-3)
#the quiesence is noe 3.6, so we just divide all of our values by this to
#  "normalize" w.r.t. 1


observed_data = (ADLEO_TWA7_c4_flux-3e-14)*10**14 +1
observed_data_err = (ADLEO_TWA7_c4_error)*10**14
popt_c4_low = (popt_c4_low)*10**14
popt_c4_high = (popt_c4_high)*10**14

# observed_data = (ADLEO_TWA7_c4_flux/3.6)*10**14 
# observed_data_err = (ADLEO_TWA7_c4_error/3.6)*10**14
# popt_c4_low = (popt_c4_low/3.6)*10**14
# popt_c4_high = (popt_c4_high/3.6)*10**14


#The following code is LEGACY and was my first attempt at caluclating these
# values. The number of 3 and 5 sigma flares each yields will be the same
#but the new method should preserve the flux/quiescence ratio while the old
# method did not.
'''
observed_data = (ADLEO_TWA7_c4_flux-np.median(ADLEO_TWA7_c4_flux))*10**14 +1
observed_data2 = (master_c4_flux-np.median(master_c4_flux))*10**14 +1
observed_data_err = (ADLEO_TWA7_c4_error)*10**14
popt_c4_low = popt_c4_low*10**14
popt_c4_high = popt_c4_high*10**14
'''


fig,subplot = plt.subplots(2, 1)

subplot[0].errorbar(master_time, observed_data, yerr = np.sqrt(observed_data_err**2+popt_c4_low[2]**2),
             ecolor = 'grey', color = 'black', elinewidth = 0.5, linewidth = 0,
             capsize = 3, marker = '.', markersize = 3, label = r'C IV (1548.2 + 1550.7)$\AA$')

subplot[0].plot(master_time,
                np.ones(len(master_time))*3*popt_c4_low[2]+1,
                color = 'orange', linewidth = 1, label = r'3$\sigma$', linestyle = 'dashed')
subplot[0].plot(master_time,
                np.ones(len(master_time))*5*popt_c4_low[2]+1,
                color = 'red', linewidth = 1, label = r'5$\sigma$')
subplot[0].plot(master_time[peaks_c4], observed_data[peaks_c4], "xr")

subplot[0].text(60,20,'low rms case', fontsize = 17, fontstyle = 'italic')
subplot[0].legend(fontsize = 13, loc = 'upper left')


subplot[1].errorbar(master_time, observed_data, yerr = np.sqrt(observed_data_err**2+popt_c4_high[2]**2),
             ecolor = 'grey', color = 'black', elinewidth = 0.5, linewidth = 0,
             capsize = 3, marker = '.', markersize = 3, label = r'C IV (1548.2 + 1550.7)$\AA$')

subplot[1].plot(master_time,
                np.ones(len(master_time))*3*popt_c4_high[2]+1,
                color = 'orange', linewidth = 1, label = r'3$\sigma$ high', linestyle = 'dashed')
subplot[1].plot(master_time,
                np.ones(len(master_time))*5*popt_c4_high[2]+1,
                color = 'red', linewidth = 1, label = r'5$\sigma$ high')
subplot[1].plot(master_time[peaks_c4], observed_data[peaks_c4], "xr")


subplot[1].text(60,20,'high rms case', fontsize = 17, fontstyle = 'italic')

fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Scaled Flux', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'TW Hya C IV Flare Propensity ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)



three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_c4)):
    temp = observed_data[peaks_c4[i]]
    if temp > 3*popt_c4_high[2]+1:
        three_sigma_list.append(i)
    if temp > 5*popt_c4_high[2]+1:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_c4[i]] + np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_high[2]**2)
    if temp > 3*popt_c4_high[2]+1:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_c4_high[2]+1:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_c4[i]] - np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_high[2]**2)
    if temp > 3*popt_c4_high[2]+1:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_c4_high[2]+1:
        five_sigma_list_err_down.append(i)

print("for C IV high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))

three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_c4)):
    temp = observed_data[peaks_c4[i]]
    if temp > 3*popt_c4_low[2]+1:
        three_sigma_list.append(i)
    if temp > 5*popt_c4_low[2]+1:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_c4[i]] + np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_low[2]**2)
    if temp > 3*popt_c4_low[2]+1:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_c4_low[2]+1:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_c4[i]] - np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_low[2]**2)
    if temp > 3*popt_c4_low[2]+1:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_c4_low[2]+1:
        five_sigma_list_err_down.append(i)

print("for C IV low:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))




#V2 of limiting case for He2

ADLEO_TWA7_He2_quies = 3.6*10**-14 #coinscidentally the same as c4
#after multiplying by 10^14 to bring the the flux up to normal magnitudes (0-3)
#the quiesence is now 3.6, so we just divide all of our values by this to
#  "normalize" w.r.t. 1

observed_data = (ADLEO_TWA7_He2_flux - 3.2e-14)*10**14+1 
observed_data_err = (ADLEO_TWA7_He2_error)*10**14
popt_He2_low = (popt_He2_low)*10**14
popt_He2_high = (popt_He2_high)*10**14

# observed_data = (ADLEO_TWA7_He2_flux/3.6)*10**14 
# observed_data_err = (ADLEO_TWA7_He2_error/3.6)*10**14
# popt_He2_low = (popt_He2_low/3.6)*10**14
# popt_He2_high = (popt_He2_high/3.6)*10**14


fig,subplot = plt.subplots(2, 1)

subplot[0].errorbar(master_time, observed_data, yerr = np.sqrt(observed_data_err**2+popt_He2_low[2]**2),
             ecolor = 'grey', color = 'black', elinewidth = 0.5, linewidth = 0,
             capsize = 3, marker = '.', markersize = 3, label = r'He II 1640$\AA$')

subplot[0].plot(master_time,
                np.ones(len(master_time))*3*popt_He2_low[2]+1,
                color = 'orange', linewidth = 1, label = r'3$\sigma$', linestyle = 'dashed')
subplot[0].plot(master_time,
                np.ones(len(master_time))*5*popt_He2_low[2]+1,
                color = 'red', linewidth = 1, label = r'5$\sigma$')
subplot[0].plot(master_time[peaks_He2], observed_data[peaks_He2], "xr")

subplot[0].text(60,6,'low rms case', fontsize = 17, fontstyle = 'italic')
subplot[0].legend(fontsize = 13, loc = 'upper left')


subplot[1].errorbar(master_time, observed_data, yerr = np.sqrt(observed_data_err**2+popt_He2_high[2]**2),
             ecolor = 'grey', color = 'black', elinewidth = 0.5, linewidth = 0,
             capsize = 3, marker = '.', markersize = 3, label = r'C IV (1548.2 + 1550.7)$\AA$')

subplot[1].plot(master_time,
                np.ones(len(master_time))*3*popt_He2_high[2]+1,
                color = 'orange', linewidth = 1, label = r'3$\sigma$ high', linestyle = 'dashed')
subplot[1].plot(master_time,
                np.ones(len(master_time))*5*popt_He2_high[2]+1,
                color = 'red', linewidth = 1, label = r'5$\sigma$ high')
subplot[1].plot(master_time[peaks_He2], observed_data[peaks_He2], "xr")


subplot[1].text(60,6,'high rms case', fontsize = 17, fontstyle = 'italic')

fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Scaled Flux', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'TW Hya He II Flare Propensity ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)



three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_He2)):
    temp = observed_data[peaks_He2[i]]
    if temp > 3*popt_He2_high[2]+1:
        three_sigma_list.append(i)
    if temp > 5*popt_He2_high[2]+1:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_He2[i]] + np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_high[2]**2)
    if temp > 3*popt_He2_high[2]+1:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_He2_high[2]+1:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_He2[i]] - np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_high[2]**2)
    if temp > 3*popt_He2_high[2]+1:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_He2_high[2]+1:
        five_sigma_list_err_down.append(i)

print("for He II high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))


three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_He2)):
    temp = observed_data[peaks_He2[i]]
    if temp > 3*popt_He2_low[2]+1:
        three_sigma_list.append(i)
    if temp > 5*popt_He2_low[2]+1:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_He2[i]] + np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_low[2]**2)
    if temp > 3*popt_He2_low[2]+1:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_He2_low[2]+1:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_He2[i]] - np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_low[2]**2)
    if temp > 3*popt_He2_low[2]+1:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_He2_low[2]+1:
        five_sigma_list_err_down.append(i)

print("for He II low:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))




# plt.figure()
# plt.plot(master_time, ADLEO_TWA7_He2_flux*10**14)


# plt.figure()
# plt.plot(master_time, ADLEO_TWA7_si4_flux*10**14/1.7)

# flax = (ADLEO_TWA7_si4_flux-np.median(ADLEO_TWA7_si4_flux))*10**14+1
# flax = (ADLEO_TWA7_si4_flux-0)*10**14+1
# plt.figure()
# plt.plot(master_time, flax)

# plt.figure()
# plt.plot(master_time, ADLEO_TWA7_si4_flux-1.9e-14+1)













'''
Now counting the number of flares with greater than 100% enhancement

'''


##############################################
#whole light curve plot with all of the species, AD Leo wtts on top of CTTS quiescence
fig,subplot = plt.subplots(3, 1)
filename_counter = 0

# for i in change_in_obs_flag:
    # flag = i in subplot0_HST_ID_fixes_list
    # if i < master_time[1117] and flag == False:
    #     subplot[0].axvline(x=i, color = 'green', alpha = 0.25, linewidth = 2)
    #     subplot[0].text(x=i,y=0.65,s=filename_list[filename_counter][0:9], color = 'green', alpha = 0.75)
    #     subplot[0].text(x=i,y=0.58,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75)
    # if i < master_time[1117]: filename_counter+=1
subplot[0].plot(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117]+1.3e-12, color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117]-ADLEO_TWA7_si4_error[0:1117]+1.3e-12,ADLEO_TWA7_si4_flux[0:1117]+ ADLEO_TWA7_si4_error[0:1117]+1.3e-12,color = 'green', alpha = 0.1)

subplot[0].plot(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]+4e-12, color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]-ADLEO_TWA7_c4_error[0:1117]+4e-12,ADLEO_TWA7_c4_flux[0:1117]+ ADLEO_TWA7_c4_error[0:1117]+4e-12,color = 'blue', alpha = 0.1)

subplot[0].plot(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117]+1.7e-12, color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117]-ADLEO_TWA7_He2_error[0:1117]+1.7e-12,ADLEO_TWA7_He2_flux[0:1117]+ ADLEO_TWA7_He2_error[0:1117]+1.7e-12,color = 'red', alpha = 0.1)




subplot[0].legend(fontsize = 13, loc = 'upper left')
# subplot[0].set_ylim(0,1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234]+1.3e-12, color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234]-ADLEO_TWA7_si4_error[1117:2234]+1.3e-12,ADLEO_TWA7_si4_flux[1117:2234]+ ADLEO_TWA7_si4_error[1117:2234]+1.3e-12,color = 'green', alpha = 0.1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]+4e-12, color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]-ADLEO_TWA7_c4_error[1117:2234]+4e-12,ADLEO_TWA7_c4_flux[1117:2234]+ ADLEO_TWA7_c4_error[1117:2234]+4e-12,color = 'blue', alpha = 0.1)

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234]+1.7e-12, color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234]-ADLEO_TWA7_He2_error[1117:2234]+1.7e-12,ADLEO_TWA7_He2_flux[1117:2234]+ ADLEO_TWA7_He2_error[1117:2234]+1.7e-12,color = 'red', alpha = 0.1)




subplot[2].plot(master_time[2234:], ADLEO_TWA7_si4_flux[2234:]+1.3e-12, color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_si4_flux[2234:]-ADLEO_TWA7_si4_error[2234:]+1.3e-12,ADLEO_TWA7_si4_flux[2234:]+ ADLEO_TWA7_si4_error[2234:]+1.3e-12,color = 'green', alpha = 0.1)

subplot[2].plot(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]+4e-12, color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]-ADLEO_TWA7_c4_error[2234:]+4e-12,ADLEO_TWA7_c4_flux[2234:]+ ADLEO_TWA7_c4_error[2234:]+4e-12,color = 'blue', alpha = 0.1)

subplot[2].plot(master_time[2234:], ADLEO_TWA7_He2_flux[2234:]+1.7e-12, color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_He2_flux[2234:]-ADLEO_TWA7_He2_error[2234:]+1.7e-12,ADLEO_TWA7_He2_flux[2234:]+ ADLEO_TWA7_He2_error[2234:]+1.7e-12,color = 'red', alpha = 0.1)




# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'$Flux(AD Leo WTTS) + F_q (CTTS)$  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to quiescent level of TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 
#highs
subplot[0].axhline(y=1.3e-12, color = 'green', alpha = 0.25, linewidth = 2)
subplot[0].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)
subplot[0].axhline(y=1.7e-12, color = 'red', alpha = 0.25, linewidth = 2)

subplot[1].axhline(y=1.3e-12, color = 'green', alpha = 0.25, linewidth = 2)
subplot[1].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)
subplot[1].axhline(y=1.7e-12, color = 'red', alpha = 0.25, linewidth = 2)

subplot[2].axhline(y=1.3e-12, color = 'green', alpha = 0.25, linewidth = 2)
subplot[2].axhline(y=4e-12, color = 'blue', alpha = 0.25, linewidth = 2)
subplot[2].axhline(y=1.7e-12, color = 'red', alpha = 0.25, linewidth = 2)













plt.close('all')
##############################################
#setting these AGAIN
#min/max limiting scenario
popt_si4_low = (1*10**(-12))*np.array([0, 0, 0.0139226519]) # computed in gaussian stats file
popt_si4_high = (1*10**(-12))*np.array([0,1.3 , 0.0444758789]) # computed in gaussian stats file

popt_c4_low = (1*10**(-12))*np.array([6.74862555, 1.535871794 , 0.048451233]) # computed in gaussian stats file
popt_c4_high = (1*10**(-12))*np.array([7.34673491, 4.0389487 , 0.06856082]) # computed in gaussian stats file

popt_He2_low = (1*10**(-12))*np.array([0,0, 0.0284057926]) # computed in gaussian stats file
popt_He2_high = (1*10**(-12))*np.array([0,1.7, 0.05246809048]) # computed in gaussian stats file

print('\n');print('\n');print('\n')
print('These are THE current estimations for the number of flares we expect to see above' 
      'the CTTS RMS when AD Leo is scaled to the level of the WTTS and and added'
       'to the CTTS signal: \n \n')






################################################
#Si IV light curve plot with all of the species, AD Leo wtts on top of CTTS quiescence
#finding delta F,q
# Fq_wtts_m = np.mean(ADLEO_TWA7_si4_flux[125:151])  #corresponding time is 2.4 to 3.0 ks
Fq_wtts_m = np.mean(ADLEO_TWA7_si4_flux[1000:1050])  #corresponding time is 20 to 21 ks
observed_data = ADLEO_TWA7_si4_flux+popt_si4_high[1] - Fq_wtts_m
observed_data_err = ADLEO_TWA7_si4_error



fig,subplot = plt.subplots(3, 1)
filename_counter = 0


subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_si4_high[1], color = 'black', alpha = .5, linewidth = 2)
# subplot[0].fill_between(master_time[0:1117], popt_si4_high[1] - popt_si4_high[2], popt_si4_high[1] + popt_si4_high[2], alpha = 0.25, color = 'black')
subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_si4_high[1]+3*popt_si4_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_si4_high[1]+5*popt_si4_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')


subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_si4_high[1], color = 'black', alpha = .5, linewidth = 2,label = '1$\sigma$')
# subplot[1].fill_between(master_time[1117:2234], popt_si4_high[1] - popt_si4_high[2], popt_si4_high[1] + popt_si4_high[2], alpha = 0.25, color = 'black')
subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_si4_high[1]+3*popt_si4_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_si4_high[1]+5*popt_si4_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')


subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_si4_high[1], color = 'black', alpha = .5, linewidth = 2,label = '1$\sigma$')
# subplot[2].fill_between(master_time[2234:], popt_si4_high[1] - popt_si4_high[2], popt_si4_high[1] + popt_si4_high[2], alpha = 0.25, color = 'black')
subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_si4_high[1]+3*popt_si4_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_si4_high[1]+5*popt_si4_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')






subplot[0].plot(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117]+popt_si4_high[1] - Fq_wtts_m, color = 'green', linewidth = 2,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_si4_flux[0:1117]-np.sqrt(ADLEO_TWA7_si4_error[0:1117]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m,ADLEO_TWA7_si4_flux[0:1117]+ np.sqrt(ADLEO_TWA7_si4_error[0:1117]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m,color = 'green', alpha = 0.25)
subplot[0].plot(master_time[peaks_si4], observed_data[peaks_si4], "xr")

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234]+popt_si4_high[1] - Fq_wtts_m, color = 'green', linewidth = 2,label = r'Si IV (1394 + 1403)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_si4_flux[1117:2234]-np.sqrt(ADLEO_TWA7_si4_error[1117:2234]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m,ADLEO_TWA7_si4_flux[1117:2234]+ np.sqrt(ADLEO_TWA7_si4_error[1117:2234]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m,color = 'green', alpha = 0.25)
subplot[1].plot(master_time[peaks_si4], observed_data[peaks_si4], "xr")


subplot[2].plot(master_time[2234:], ADLEO_TWA7_si4_flux[2234:]+popt_si4_high[1] - Fq_wtts_m, color = 'green', linewidth = 2,label = r'Si IV (1394 + 1403)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_si4_flux[2234:]-np.sqrt(ADLEO_TWA7_si4_error[2234:]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m,ADLEO_TWA7_si4_flux[2234:]+ np.sqrt(ADLEO_TWA7_si4_error[2234:]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m,color = 'green', alpha = 0.25)
subplot[2].plot(master_time[peaks_si4], observed_data[peaks_si4], "xr")


# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'$F_{WTTS, \/\/M^*} - F_{q, \/\/  WTTS, \/\/  M^*} + F_{q, CTTS}$  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to quiescent level of TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 

subplot[0].set_xlim(0,22.32)
subplot[0].set_ylim(1.25e-12,1.65e-12)
subplot[1].set_xlim(22.34,44.66)
subplot[1].set_ylim(1.25e-12,1.65e-12)
subplot[2].set_xlim(44.66,66.98)
subplot[2].set_ylim(1.25e-12,1.65e-12)
subplot[1].set_ylim(1.25e-12,1.65e-12)
subplot[2].set_ylim(1.25e-12,1.65e-12)
subplot[0].legend(fontsize = 13, loc = 'upper left')





three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_si4)):
    temp = observed_data[peaks_si4[i]]
    if temp > 3*popt_si4_high[2]+popt_si4_high[1]:
        three_sigma_list.append(i)
    if temp > 5*popt_si4_high[2]+popt_si4_high[1]:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_si4[i]] + np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_high[2]**2)
    if temp > 3*popt_si4_high[2]+popt_si4_high[1]:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_si4_high[2]+popt_si4_high[1]:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_si4[i]] - np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_high[2]**2)
    if temp > 3*popt_si4_high[2]+popt_si4_high[1]:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_si4_high[2]+popt_si4_high[1]:
        five_sigma_list_err_down.append(i)

print("for Si IV high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))



################################################
#C IV light curve plot, AD Leo wtts on top of CTTS quiescence
#finding delta F,q
# Fq_wtts_m = np.mean(ADLEO_TWA7_c4_flux[125:151])  #corresponding time is 2.4 to 3.0 ks
Fq_wtts_m = np.mean(ADLEO_TWA7_c4_flux[1000:1050])  #corresponding time is 20 to 21 ks
observed_data = ADLEO_TWA7_c4_flux+popt_c4_high[1] - Fq_wtts_m
observed_data_err = ADLEO_TWA7_c4_error


fig,subplot = plt.subplots(3, 1)
filename_counter = 0


subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_c4_high[1], color = 'black', alpha = .5, linewidth = 2)
# subplot[0].fill_between(master_time[0:1117], popt_c4_high[1] - popt_c4_high[2], popt_c4_high[1] + popt_c4_high[2], alpha = 0.25, color = 'black')
subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_c4_high[1]+3*popt_c4_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_c4_high[1]+5*popt_c4_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')


subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_c4_high[1], color = 'black', alpha = .5, linewidth = 2,label = '1$\sigma$')
# subplot[1].fill_between(master_time[1117:2234], popt_c4_high[1] - popt_c4_high[2], popt_c4_high[1] + popt_c4_high[2], alpha = 0.25, color = 'black')
subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_c4_high[1]+3*popt_c4_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_c4_high[1]+5*popt_c4_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')


subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_c4_high[1], color = 'black', alpha = .5, linewidth = 2,label = '1$\sigma$')
# subplot[2].fill_between(master_time[2234:], popt_c4_high[1] - popt_c4_high[2], popt_c4_high[1] + popt_c4_high[2], alpha = 0.25, color = 'black')
subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_c4_high[1]+3*popt_c4_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_c4_high[1]+5*popt_c4_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')






subplot[0].plot(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]+popt_c4_high[1] - Fq_wtts_m, color = 'blue', linewidth = 2,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_c4_flux[0:1117]-np.sqrt(ADLEO_TWA7_c4_error[0:1117]**2+popt_c4_high[2]**2)+popt_c4_high[1] - Fq_wtts_m,ADLEO_TWA7_c4_flux[0:1117]+ np.sqrt(ADLEO_TWA7_c4_error[0:1117]**2+popt_c4_high[2]**2)+popt_c4_high[1] - Fq_wtts_m,color = 'blue', alpha = 0.25)
subplot[0].plot(master_time[peaks_c4], observed_data[peaks_c4], "xr")

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]+popt_c4_high[1] - Fq_wtts_m, color = 'blue', linewidth = 2,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_c4_flux[1117:2234]-np.sqrt(ADLEO_TWA7_c4_error[1117:2234]**2+popt_c4_high[2]**2)+popt_c4_high[1] - Fq_wtts_m,ADLEO_TWA7_c4_flux[1117:2234]+ np.sqrt(ADLEO_TWA7_c4_error[1117:2234]**2+popt_c4_high[2]**2)+popt_c4_high[1] - Fq_wtts_m,color = 'blue', alpha = 0.25)
subplot[1].plot(master_time[peaks_c4], observed_data[peaks_c4], "xr")


subplot[2].plot(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]+popt_c4_high[1] - Fq_wtts_m, color = 'blue', linewidth = 2,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_c4_flux[2234:]-np.sqrt(ADLEO_TWA7_c4_error[2234:]**2+popt_c4_high[2]**2)+popt_c4_high[1] - Fq_wtts_m,ADLEO_TWA7_c4_flux[2234:]+ np.sqrt(ADLEO_TWA7_c4_error[2234:]**2+popt_c4_high[2]**2)+popt_c4_high[1] - Fq_wtts_m,color = 'blue', alpha = 0.25)
subplot[2].plot(master_time[peaks_c4], observed_data[peaks_c4], "xr")


# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'$F_{WTTS, \/\/M^*} - F_{q, \/\/  WTTS, \/\/  M^*} + F_{q, CTTS}$  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to quiescent level of TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 

subplot[0].set_xlim(0,22.32)
subplot[0].set_ylim(3.95e-12,4.55e-12)
subplot[1].set_xlim(22.34,44.66)
subplot[1].set_ylim(3.95e-12,4.55e-12)
subplot[2].set_xlim(44.66,66.98)
subplot[2].set_ylim(3.95e-12,4.55e-12)
subplot[1].set_ylim(3.95e-12,4.55e-12)
subplot[2].set_ylim(3.95e-12,4.55e-12)
subplot[0].legend(fontsize = 13, loc = 'upper left')


three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_c4)):
    temp = observed_data[peaks_c4[i]]
    if temp > 3*popt_c4_high[2]+popt_c4_high[1]:
        three_sigma_list.append(i)
    if temp > 5*popt_c4_high[2]+popt_c4_high[1]:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_c4[i]] + np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_high[2]**2)
    if temp > 3*popt_c4_high[2]+popt_c4_high[1]:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_c4_high[2]+popt_c4_high[1]:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_c4[i]] - np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_high[2]**2)
    if temp > 3*popt_c4_high[2]+popt_c4_high[1]:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_c4_high[2]+popt_c4_high[1]:
        five_sigma_list_err_down.append(i)

print("for C IV high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))





################################################
#He II light curve plot with all of the species, AD Leo wtts on top of CTTS quiescence
#finding delta F,q
# Fq_wtts_m = np.mean(ADLEO_TWA7_He2_flux[125:151])  #corresponding time is 2.4 to 3.0 ks
Fq_wtts_m = np.mean(ADLEO_TWA7_He2_flux[1000:1050])  #corresponding time is 20 to 21 ks

observed_data = ADLEO_TWA7_He2_flux+popt_He2_high[1] - Fq_wtts_m
observed_data_err = ADLEO_TWA7_He2_error


fig,subplot = plt.subplots(3, 1)
filename_counter = 0


subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_He2_high[1], color = 'black', alpha = .5, linewidth = 2)
# subplot[0].fill_between(master_time[0:1117], popt_He2_high[1] - popt_He2_high[2], popt_He2_high[1] + popt_He2_high[2], alpha = 0.25, color = 'black')
subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_He2_high[1]+3*popt_He2_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
# subplot[0].plot(master_time[0:1117],np.ones(len(master_time[0:1117]))*popt_He2_high[1]+5*popt_He2_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')


subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_He2_high[1], color = 'black', alpha = .5, linewidth = 2,label = '1$\sigma$')
# subplot[1].fill_between(master_time[1117:2234], popt_He2_high[1] - popt_He2_high[2], popt_He2_high[1] + popt_He2_high[2], alpha = 0.25, color = 'black')
subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_He2_high[1]+3*popt_He2_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[1].plot(master_time[1117:2234],np.ones(len(master_time[1117:2234]))*popt_He2_high[1]+5*popt_He2_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')


subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_He2_high[1], color = 'black', alpha = .5, linewidth = 2,label = '1$\sigma$')
# subplot[2].fill_between(master_time[2234:], popt_He2_high[1] - popt_He2_high[2], popt_He2_high[1] + popt_He2_high[2], alpha = 0.25, color = 'black')
subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_He2_high[1]+3*popt_He2_high[2], color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[2].plot(master_time[2234:],np.ones(len(master_time[2234:]))*popt_He2_high[1]+5*popt_He2_high[2], color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')






subplot[0].plot(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117]+popt_He2_high[1] - Fq_wtts_m, color = 'black',alpha = 0.5, linewidth = 2,label = r'He II 1640$\AA$')
subplot[0].fill_between(master_time[0:1117], ADLEO_TWA7_He2_flux[0:1117]-np.sqrt(ADLEO_TWA7_He2_error[0:1117]**2+popt_He2_high[2]**2)+popt_He2_high[1] - Fq_wtts_m,ADLEO_TWA7_He2_flux[0:1117]+ np.sqrt(ADLEO_TWA7_He2_error[0:1117]**2+popt_He2_high[2]**2)+popt_He2_high[1] - Fq_wtts_m,color = 'black', alpha = 0.25)
subplot[0].plot(master_time[peaks_He2], observed_data[peaks_He2], "xr")

subplot[1].plot(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234]+popt_He2_high[1] - Fq_wtts_m, color = 'black',alpha = 0.5, linewidth = 2,label = r'He II 1640$\AA$')
subplot[1].fill_between(master_time[1117:2234], ADLEO_TWA7_He2_flux[1117:2234]-np.sqrt(ADLEO_TWA7_He2_error[1117:2234]**2+popt_He2_high[2]**2)+popt_He2_high[1] - Fq_wtts_m,ADLEO_TWA7_He2_flux[1117:2234]+ np.sqrt(ADLEO_TWA7_He2_error[1117:2234]**2+popt_He2_high[2]**2)+popt_He2_high[1] - Fq_wtts_m,color = 'black', alpha = 0.25)
subplot[1].plot(master_time[peaks_He2], observed_data[peaks_He2], "xr")


subplot[2].plot(master_time[2234:], ADLEO_TWA7_He2_flux[2234:]+popt_He2_high[1] - Fq_wtts_m, color = 'black',alpha = 0.5, linewidth = 2,label = r'He II 1640$\AA$')
subplot[2].fill_between(master_time[2234:], ADLEO_TWA7_He2_flux[2234:]-np.sqrt(ADLEO_TWA7_He2_error[2234:]**2+popt_He2_high[2]**2)+popt_He2_high[1] - Fq_wtts_m,ADLEO_TWA7_He2_flux[2234:]+ np.sqrt(ADLEO_TWA7_He2_error[2234:]**2+popt_He2_high[2]**2)+popt_He2_high[1] - Fq_wtts_m,color = 'black', alpha = 0.25)
subplot[2].plot(master_time[peaks_He2], observed_data[peaks_He2], "xr")


# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'$F_{WTTS, \/\/M^*} - F_{q, \/\/  WTTS, \/\/  M^*} + F_{q, CTTS}$  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to quiescent level of TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 

subplot[0].set_xlim(0,22.32)
subplot[0].set_ylim(1.6e-12,1.9e-12)
subplot[1].set_xlim(22.34,44.66)
subplot[1].set_ylim(1.6e-12,1.9e-12)
subplot[2].set_xlim(44.66,66.98)
subplot[2].set_ylim(1.6e-12,1.9e-12)
subplot[1].set_ylim(1.6e-12,1.9e-12)
subplot[2].set_ylim(1.6e-12,1.9e-12)
subplot[0].legend(fontsize = 13, loc = 'upper left')


three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
for i in range(len(peaks_He2)):
    temp = observed_data[peaks_He2[i]]
    if temp > 3*popt_He2_high[2]+popt_He2_high[1]:
        three_sigma_list.append(i)
    if temp > 5*popt_He2_high[2]+popt_He2_high[1]:
        five_sigma_list.append(i)   
    temp = observed_data[peaks_He2[i]] + np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_high[2]**2)
    if temp > 3*popt_He2_high[2]+popt_He2_high[1]:
        three_sigma_list_err_up.append(i)
    if temp > 5*popt_He2_high[2]+popt_He2_high[1]:
        five_sigma_list_err_up.append(i)
    temp = observed_data[peaks_He2[i]] - np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_high[2]**2)
    if temp > 3*popt_He2_high[2]+popt_He2_high[1]:
        three_sigma_list_err_down.append(i)
    if temp > 5*popt_He2_high[2]+popt_He2_high[1]:
        five_sigma_list_err_down.append(i)

print("for He II high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))







######
print("now computing the estimates using the Hilton 2011 approach:")
three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
Fq_wtts_m = np.mean(ADLEO_TWA7_si4_flux[1000:1050])  #corresponding time is 20 to 21 ks
observed_data = ADLEO_TWA7_si4_flux+popt_si4_high[1] - Fq_wtts_m
observed_data_err = ADLEO_TWA7_si4_error
for i in range(len(peaks_si4)):
    #finding base 3 sigma events
    temp = observed_data[peaks_si4[i]]
    if temp > 3*popt_si4_high[2]+popt_si4_high[1]:
        temp2 = observed_data[peaks_si4[i]+1];temp3 = observed_data[peaks_si4[i]+2];
        temp4 = observed_data[peaks_si4[i]-1];temp5 = observed_data[peaks_si4[i]-2];
        if temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp3 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list.append(i)
        elif temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list.append(i)
        elif temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp5 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list.append(i)
    #finding base 5 sigma events
    if temp > 5*popt_si4_high[2]+popt_si4_high[1]:
        temp2 = observed_data[peaks_si4[i]+1];temp3 = observed_data[peaks_si4[i]+2];
        temp4 = observed_data[peaks_si4[i]-1];temp5 = observed_data[peaks_si4[i]-2];
        if temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp3 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list.append(i) 
        elif temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list.append(i)
        elif temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp5 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list.append(i)
    #finding upper 3 sigma events
    temp_err = np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_high[2]**2)
    temp = observed_data[peaks_si4[i]] + temp_err
    if temp > 3*popt_si4_high[2]+popt_si4_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_si4[i]+1]**2 + popt_si4_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_si4[i]+2]**2 + popt_si4_high[2]**2)
        temp_err4 = np.sqrt(observed_data_err[peaks_si4[i]-1]**2 + popt_si4_high[2]**2)
        temp_err5 = np.sqrt(observed_data_err[peaks_si4[i]-2]**2 + popt_si4_high[2]**2)        
        temp2 = observed_data[peaks_si4[i]+1] + temp_err2;
        temp3 = observed_data[peaks_si4[i]+2] + temp_err3;
        temp4 = observed_data[peaks_si4[i]-1] + temp_err4;
        temp5 = observed_data[peaks_si4[i]-2] + temp_err5;
        if temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp3 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list_err_up.append(i)
        elif temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list_err_up.append(i)
        elif temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp5 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list_err_up.append(i)
        
    #finding lower 3 sigma events            
    temp_err = np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_high[2]**2)
    temp = observed_data[peaks_si4[i]] - temp_err
    if temp > 3*popt_si4_high[2]+popt_si4_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_si4[i]+1]**2 + popt_si4_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_si4[i]+2]**2 + popt_si4_high[2]**2)
        temp_err4 = np.sqrt(observed_data_err[peaks_si4[i]-1]**2 + popt_si4_high[2]**2)
        temp_err5 = np.sqrt(observed_data_err[peaks_si4[i]-2]**2 + popt_si4_high[2]**2)        
        temp2 = observed_data[peaks_si4[i]+1] - temp_err2;
        temp3 = observed_data[peaks_si4[i]+2] - temp_err3;
        temp4 = observed_data[peaks_si4[i]-1] - temp_err4;
        temp5 = observed_data[peaks_si4[i]-2] - temp_err5;
        if temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp3 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list_err_down.append(i)
        elif temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list_err_down.append(i)
        elif temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp5 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            three_sigma_list_err_down.append(i)
    #finding upper 5 sigma events
    temp_err = np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_high[2]**2)
    temp = observed_data[peaks_si4[i]] + temp_err
    if temp > 5*popt_si4_high[2]+popt_si4_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_si4[i]+1]**2 + popt_si4_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_si4[i]+2]**2 + popt_si4_high[2]**2)
        temp_err4 = np.sqrt(observed_data_err[peaks_si4[i]-1]**2 + popt_si4_high[2]**2)
        temp_err5 = np.sqrt(observed_data_err[peaks_si4[i]-2]**2 + popt_si4_high[2]**2)         
        temp2 = observed_data[peaks_si4[i]+1] + temp_err2;
        temp3 = observed_data[peaks_si4[i]+2] + temp_err3;
        temp4 = observed_data[peaks_si4[i]-1] + temp_err4;
        temp5 = observed_data[peaks_si4[i]-2] + temp_err5;
        if temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp3 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list_err_up.append(i)
        elif temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list_err_up.append(i)
        elif temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp5 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list_err_up.append(i)
    #finding lower 5 sigma events        
    temp_err = np.sqrt(observed_data_err[peaks_si4[i]]**2 + popt_si4_high[2]**2)
    temp = observed_data[peaks_si4[i]] - temp_err
    if temp > 5*popt_si4_high[2]+popt_si4_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_si4[i]+1]**2 + popt_si4_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_si4[i]+2]**2 + popt_si4_high[2]**2)
        temp_err4 = np.sqrt(observed_data_err[peaks_si4[i]-1]**2 + popt_si4_high[2]**2)
        temp_err5 = np.sqrt(observed_data_err[peaks_si4[i]-2]**2 + popt_si4_high[2]**2)         
        temp2 = observed_data[peaks_si4[i]+1] - temp_err2;
        temp3 = observed_data[peaks_si4[i]+2] - temp_err3;
        temp4 = observed_data[peaks_si4[i]-1] - temp_err4;
        temp5 = observed_data[peaks_si4[i]-2] - temp_err5;
        if temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp3 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list_err_down.append(i)
        elif temp2 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list_err_down.append(i)
        elif temp4 > 2.5*popt_si4_high[2]+popt_si4_high[1] and temp5 > 2.5*popt_si4_high[2]+popt_si4_high[1]:
            five_sigma_list_err_down.append(i)
print("for S IV high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))




three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
Fq_wtts_m = np.mean(ADLEO_TWA7_c4_flux[1000:1050])  #corresponding time is 20 to 21 ks
observed_data = ADLEO_TWA7_c4_flux+popt_c4_high[1] - Fq_wtts_m
for i in range(len(peaks_c4)):
    #finding base 3 sigma events
    temp = observed_data[peaks_c4[i]]
    if temp > 3*popt_c4_high[2]+popt_c4_high[1]:
        temp2 = observed_data[peaks_c4[i]+1];temp3 = observed_data[peaks_c4[i]+2];
        if temp2 > 2.5*popt_c4_high[2]+popt_c4_high[1] and temp3 > 2.5*popt_c4_high[2]+popt_c4_high[1]:
            three_sigma_list.append(i)
    #finding base 5 sigma events
    if temp > 5*popt_c4_high[2]+popt_c4_high[1]:
        temp2 = observed_data[peaks_c4[i]+1];temp3 = observed_data[peaks_c4[i]+2];
        if temp2 > 2.5*popt_c4_high[2]+popt_c4_high[1] and temp3 > 2.5*popt_c4_high[2]+popt_c4_high[1]:
            five_sigma_list.append(i)    
    #finding upper 3 sigma events
    temp_err = np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_high[2]**2)
    temp = observed_data[peaks_c4[i]] + temp_err
    if temp > 3*popt_c4_high[2]+popt_c4_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_c4[i]+1]**2 + popt_c4_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_c4[i]+2]**2 + popt_c4_high[2]**2)        
        temp2 = observed_data[peaks_c4[i]+1] + temp_err2;
        temp3 = observed_data[peaks_c4[i]+2] + temp_err3;
        if temp2 > 2.5*popt_c4_high[2]+popt_c4_high[1] and temp3 > 2.5*popt_c4_high[2]+popt_c4_high[1]:
            three_sigma_list_err_up.append(i)
    #finding lower 3 sigma events            
    temp_err = np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_high[2]**2)
    temp = observed_data[peaks_c4[i]] - temp_err
    if temp > 3*popt_c4_high[2]+popt_c4_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_c4[i]+1]**2 + popt_c4_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_c4[i]+2]**2 + popt_c4_high[2]**2)        
        temp2 = observed_data[peaks_c4[i]+1] - temp_err2;
        temp3 = observed_data[peaks_c4[i]+2] - temp_err3;
        if temp2 > 2.5*popt_c4_high[2]+popt_c4_high[1] and temp3 > 2.5*popt_c4_high[2]+popt_c4_high[1]:
            three_sigma_list_err_down.append(i)
    #finding upper 5 sigma events
    temp_err = np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_high[2]**2)
    temp = observed_data[peaks_c4[i]] + temp_err
    if temp > 5*popt_c4_high[2]+popt_c4_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_c4[i]+1]**2 + popt_c4_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_c4[i]+2]**2 + popt_c4_high[2]**2)        
        temp2 = observed_data[peaks_c4[i]+1] + temp_err2;
        temp3 = observed_data[peaks_c4[i]+2] + temp_err3;
        if temp2 > 2.5*popt_c4_high[2]+popt_c4_high[1] and temp3 > 2.5*popt_c4_high[2]+popt_c4_high[1]:
            five_sigma_list_err_up.append(i)
    #finding lower 5 sigma events        
    temp_err = np.sqrt(observed_data_err[peaks_c4[i]]**2 + popt_c4_high[2]**2)
    temp = observed_data[peaks_c4[i]] - temp_err
    if temp > 5*popt_c4_high[2]+popt_c4_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_c4[i]+1]**2 + popt_c4_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_c4[i]+2]**2 + popt_c4_high[2]**2)        
        temp2 = observed_data[peaks_c4[i]+1] - temp_err2;
        temp3 = observed_data[peaks_c4[i]+2] - temp_err3;
        if temp2 > 2.5*popt_c4_high[2]+popt_c4_high[1] and temp3 > 2.5*popt_c4_high[2]+popt_c4_high[1]:
            five_sigma_list_err_down.append(i)
print("for C IV high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))




three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
Fq_wtts_m = np.mean(ADLEO_TWA7_He2_flux[1000:1050])  #corresponding time is 20 to 21 ks
observed_data = ADLEO_TWA7_He2_flux+popt_He2_high[1] - Fq_wtts_m
for i in range(len(peaks_He2)):
    #finding base 3 sigma events
    temp = observed_data[peaks_He2[i]]
    if temp > 3*popt_He2_high[2]+popt_He2_high[1]:
        temp2 = observed_data[peaks_He2[i]+1];temp3 = observed_data[peaks_He2[i]+2];
        if temp2 > 2.5*popt_He2_high[2]+popt_He2_high[1] and temp3 > 2.5*popt_He2_high[2]+popt_He2_high[1]:
            three_sigma_list.append(i)
    #finding base 5 sigma events
    if temp > 5*popt_He2_high[2]+popt_He2_high[1]:
        temp2 = observed_data[peaks_He2[i]+1];temp3 = observed_data[peaks_He2[i]+2];
        if temp2 > 2.5*popt_He2_high[2]+popt_He2_high[1] and temp3 > 2.5*popt_He2_high[2]+popt_He2_high[1]:
            five_sigma_list.append(i)    
    #finding upper 3 sigma events
    temp_err = np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_high[2]**2)
    temp = observed_data[peaks_He2[i]] + temp_err
    if temp > 3*popt_He2_high[2]+popt_He2_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_He2[i]+1]**2 + popt_He2_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_He2[i]+2]**2 + popt_He2_high[2]**2)        
        temp2 = observed_data[peaks_He2[i]+1] + temp_err2;
        temp3 = observed_data[peaks_He2[i]+2] + temp_err3;
        if temp2 > 2.5*popt_He2_high[2]+popt_He2_high[1] and temp3 > 2.5*popt_He2_high[2]+popt_He2_high[1]:
            three_sigma_list_err_up.append(i)
    #finding lower 3 sigma events            
    temp_err = np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_high[2]**2)
    temp = observed_data[peaks_He2[i]] - temp_err
    if temp > 3*popt_He2_high[2]+popt_He2_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_He2[i]+1]**2 + popt_He2_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_He2[i]+2]**2 + popt_He2_high[2]**2)        
        temp2 = observed_data[peaks_He2[i]+1] - temp_err2;
        temp3 = observed_data[peaks_He2[i]+2] - temp_err3;
        if temp2 > 2.5*popt_He2_high[2]+popt_He2_high[1] and temp3 > 2.5*popt_He2_high[2]+popt_He2_high[1]:
            three_sigma_list_err_down.append(i)
    #finding upper 5 sigma events
    temp_err = np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_high[2]**2)
    temp = observed_data[peaks_He2[i]] + temp_err
    if temp > 5*popt_He2_high[2]+popt_He2_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_He2[i]+1]**2 + popt_He2_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_He2[i]+2]**2 + popt_He2_high[2]**2)        
        temp2 = observed_data[peaks_He2[i]+1] + temp_err2;
        temp3 = observed_data[peaks_He2[i]+2] + temp_err3;
        if temp2 > 2.5*popt_He2_high[2]+popt_He2_high[1] and temp3 > 2.5*popt_He2_high[2]+popt_He2_high[1]:
            five_sigma_list_err_up.append(i)
    #finding lower 5 sigma events        
    temp_err = np.sqrt(observed_data_err[peaks_He2[i]]**2 + popt_He2_high[2]**2)
    temp = observed_data[peaks_He2[i]] - temp_err
    if temp > 5*popt_He2_high[2]+popt_He2_high[1]:
        temp_err2 = np.sqrt(observed_data_err[peaks_He2[i]+1]**2 + popt_He2_high[2]**2)
        temp_err3 = np.sqrt(observed_data_err[peaks_He2[i]+2]**2 + popt_He2_high[2]**2)        
        temp2 = observed_data[peaks_He2[i]+1] - temp_err2;
        temp3 = observed_data[peaks_He2[i]+2] - temp_err3;
        if temp2 > 2.5*popt_He2_high[2]+popt_He2_high[1] and temp3 > 2.5*popt_He2_high[2]+popt_He2_high[1]:
            five_sigma_list_err_down.append(i)
print("for He II high:")
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))





#################   c4 accretion v flare energy contributions ##############

total_time =master_time[-1].copy()*10**3#convert ks to s
distance_correction = 4*np.pi*(56*3.086e18)**2 #56 p to cm




# Fq_wtts_m_c4 = np.mean(ADLEO_TWA7_c4_flux[125:151]);
# std_Fq_wtts_m_c4 = np.std(ADLEO_TWA7_c4_flux[125:151])
Fq_wtts_m_c4 = np.mean(ADLEO_TWA7_c4_flux[1000:1050]); #20 to 21ks
std_Fq_wtts_m_c4 = np.std(ADLEO_TWA7_c4_flux[1000:1050])
accretion_contribution_list = []
flare_contribution_list = []
for i in tqdm(range(10000)):
    flare_signal_c4 = []
    CTTS_signal = np.ones(len(master_time))*np.random.normal(popt_c4_high[1], popt_c4_high[2])
    WTTS_signal = np.ones(len(master_time))*np.random.normal(Fq_wtts_m_c4, std_Fq_wtts_m_c4)
    for j in range(len(ADLEO_TWA7_c4_flux)):
        temp_mu = ADLEO_TWA7_c4_flux[j];temp_std = ADLEO_TWA7_c4_error[j];
        temp_flare_signal = np.random.normal(temp_mu, temp_std)-np.random.normal(Fq_wtts_m_c4, std_Fq_wtts_m_c4)
        flare_signal_c4.append(temp_flare_signal)

    total_CTTS_energy_c4 = np.trapz(CTTS_signal,master_time)*10**3*distance_correction #in ergs
    WTTS_q_energy_c4     = np.trapz(WTTS_signal,master_time)*10**3*distance_correction #in ergs
    flare_signal_energy_c4 = np.trapz(flare_signal_c4,master_time)*10**3*distance_correction
    accretion_energy_c4 = total_CTTS_energy_c4 - WTTS_q_energy_c4
    
    # flare_signal_c4 = (ADLEO_TWA7_c4_flux - Fq_wtts_m) 
    # total_CTTS_energy_c4 = popt_c4_high[1]*total_time*distance_correction #in ergs
    # WTTS_q_energy_c4 = Fq_wtts_m*total_time*distance_correction #in ergs
    # flare_signal_energy_c4 = np.trapz(flare_signal_c4,master_time)*10**3*distance_correction
    # accretion_energy_c4 = total_CTTS_energy_c4 - WTTS_q_energy_c4 - flare_signal_energy_c4
    accretion_contribution_list.append(accretion_energy_c4/total_CTTS_energy_c4)
    flare_contribution_list.append(flare_signal_energy_c4/total_CTTS_energy_c4)

print('\n');print('\n');print('\n')
print('Flare contribution = %s +/- %s' % (np.mean(flare_contribution_list)*100,100*np.std(flare_contribution_list)))
print('Accretion contribution = %s +/- %s' % (np.mean(accretion_contribution_list)*100,100*np.std(accretion_contribution_list)))


# print('\n');print('\n');print('\n')
# temp = str(flare_signal_energy_c4)
# print('Total energy in CIV flares: %s ergs' % temp)
# temp =str(accretion_energy_c4)
# print('Total energy in CIV attributed to accretion: %s ergs' % temp)
# temp = str(flare_signal_energy_c4/total_CTTS_energy_c4) 
# print('Fraction of energy in C IV flares: %s' % temp[0:6])
# temp = str(accretion_energy_c4/total_CTTS_energy_c4)
# print('Fraction of energy in C IV attributed to accretion: %s' % temp[0:6])

##############################################################################


plt.figure()
count, bins, ignored = plt.hist(accretion_contribution_list,50, density=True)
plt.title('TW Hya Accretion Contribution Distribution')
plt.ylabel('Frequency')
plt.xlabel('Fraction Contributed')

plt.figure()
count, bins, ignored = plt.hist(flare_contribution_list,50, density=True)
plt.title('TW Hya Flare Contribution Distribution')
plt.ylabel('Frequency')
plt.xlabel('Fraction Contributed')

# plt.figure()
# plt.plot(master_time, flare_signal_c4)



#### testing stiching 67ks frames together to widen statistics #########
# Fq_wtts_m_si4 = np.mean(ADLEO_TWA7_si4_flux[1000:1050]); #20 to 21ks
# std_Fq_wtts_m_si4 = np.std(ADLEO_TWA7_si4_flux[1000:1050])
# five_sigma_list = []
# three_sigma_list = []
# num_of_67ks_frames = 1000
# for j in tqdm(range(1000)):
#     five_sigma_counter = 0
#     for k in range(num_of_67ks_frames):
#         for i in range(len(peaks_si4)):
#             potential_flare = (np.random.normal(ADLEO_TWA7_si4_flux[peaks_si4[i]], ADLEO_TWA7_si4_error[peaks_si4[i]]) -
#                                 np.random.normal(Fq_wtts_m_si4, std_Fq_wtts_m_si4)+
#                                 np.random.normal(popt_si4_high[1],popt_si4_high[2]))
#             if potential_flare > 5*popt_si4_high[2]+popt_si4_high[1]:
#                     five_sigma_counter += 1
    
#     five_sigma_list.append(five_sigma_counter)  
# plt.figure()
# count, bins, ignored = plt.hist(five_sigma_list,25, density=True)
# plt.title('TW Hya Flare Frequency')
# plt.ylabel('Frequency')
# plt.xlabel('# of flares')
# print('for S IV:')
# print(np.mean(five_sigma_list)/num_of_67ks_frames)
# print(np.std(five_sigma_list)/num_of_67ks_frames)







popt_si4_high = (1*10**(-12))*np.array([0,1.3 , 0.0444758789]) # computed in gaussian stats file
#plot for paper

Fq_wtts_m = np.mean(ADLEO_TWA7_si4_flux[1000:1050])  #corresponding time is 20 to 21 ks
observed_data = ADLEO_TWA7_si4_flux+popt_si4_high[1] - Fq_wtts_m
observed_data_err = ADLEO_TWA7_si4_error



fig,subplot = plt.subplots(3, 1)
filename_counter = 0


subplot[0].plot(master_time[0:1117],1e12*(np.ones(len(master_time[0:1117]))*popt_si4_high[1]), color = 'black', alpha = .5, linewidth = 2)
# subplot[0].fill_between(master_time[0:1117], popt_si4_high[1] - popt_si4_high[2], popt_si4_high[1] + popt_si4_high[2], alpha = 0.25, color = 'black')
subplot[0].plot(master_time[0:1117],1e12*(np.ones(len(master_time[0:1117]))*popt_si4_high[1]+3*popt_si4_high[2]), color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[0].plot(master_time[0:1117],1e12*(np.ones(len(master_time[0:1117]))*popt_si4_high[1]+5*popt_si4_high[2]), color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')


subplot[1].plot(master_time[1117:2234],1e12*(np.ones(len(master_time[1117:2234]))*popt_si4_high[1]), color = 'black', alpha = .5, linewidth = 2,label = '1$\sigma$')
# subplot[1].fill_between(master_time[1117:2234], popt_si4_high[1] - popt_si4_high[2], popt_si4_high[1] + popt_si4_high[2], alpha = 0.25, color = 'black')
subplot[1].plot(master_time[1117:2234],1e12*(np.ones(len(master_time[1117:2234]))*popt_si4_high[1]+3*popt_si4_high[2]), color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[1].plot(master_time[1117:2234],1e12*(np.ones(len(master_time[1117:2234]))*popt_si4_high[1]+5*popt_si4_high[2]), color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')


subplot[2].plot(master_time[2234:],1e12*(np.ones(len(master_time[2234:]))*popt_si4_high[1]), color = 'black', alpha = .5, linewidth = 2,label = '1$\sigma$')
# subplot[2].fill_between(master_time[2234:], popt_si4_high[1] - popt_si4_high[2], popt_si4_high[1] + popt_si4_high[2], alpha = 0.25, color = 'black')
subplot[2].plot(master_time[2234:],1e12*(np.ones(len(master_time[2234:]))*popt_si4_high[1]+3*popt_si4_high[2]), color = 'orange', alpha = 1, linewidth = 2, label = '3$\sigma$')
subplot[2].plot(master_time[2234:],1e12*(np.ones(len(master_time[2234:]))*popt_si4_high[1]+5*popt_si4_high[2]), color = 'red', alpha = 1, linewidth = 2, label = '5$\sigma$')






subplot[0].plot(master_time[0:1117], 1e12*(ADLEO_TWA7_si4_flux[0:1117]+popt_si4_high[1] - Fq_wtts_m), color = 'green', linewidth = 2,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], 1e12*(ADLEO_TWA7_si4_flux[0:1117]-np.sqrt(ADLEO_TWA7_si4_error[0:1117]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m),1e12*(ADLEO_TWA7_si4_flux[0:1117]+ np.sqrt(ADLEO_TWA7_si4_error[0:1117]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m),color = 'green', alpha = 0.25)
subplot[0].plot(master_time[peaks_si4], 1e12*observed_data[peaks_si4], "xr")

subplot[1].plot(master_time[1117:2234], 1e12*(ADLEO_TWA7_si4_flux[1117:2234]+popt_si4_high[1] - Fq_wtts_m), color = 'green', linewidth = 2,label = r'Si IV (1394 + 1403)$\AA$')
subplot[1].fill_between(master_time[1117:2234], 1e12*(ADLEO_TWA7_si4_flux[1117:2234]-np.sqrt(ADLEO_TWA7_si4_error[1117:2234]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m),1e12*(ADLEO_TWA7_si4_flux[1117:2234]+ np.sqrt(ADLEO_TWA7_si4_error[1117:2234]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m),color = 'green', alpha = 0.25)
subplot[1].plot(master_time[peaks_si4], 1e12*observed_data[peaks_si4], "xr")


subplot[2].plot(master_time[2234:], 1e12*(ADLEO_TWA7_si4_flux[2234:]+popt_si4_high[1] - Fq_wtts_m), color = 'green', linewidth = 2,label = r'Si IV (1394 + 1403)$\AA$')
subplot[2].fill_between(master_time[2234:], 1e12*(ADLEO_TWA7_si4_flux[2234:]-np.sqrt(ADLEO_TWA7_si4_error[2234:]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m),1e12*(ADLEO_TWA7_si4_flux[2234:]+ np.sqrt(ADLEO_TWA7_si4_error[2234:]**2+popt_si4_high[2]**2)+popt_si4_high[1] - Fq_wtts_m),color = 'green', alpha = 0.25)
subplot[2].plot(master_time[peaks_si4], 1e12*observed_data[peaks_si4], "xr")


# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'$F_{WTTS, \/\/M^*} - F_{q, \/\/  WTTS, \/\/  M^*} + F_{q, CTTS}$  ($erg \ cm^{-2} s^{-1}*10^{-12}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'AD Leo Scaled to quiescent level of TWA-7 ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 

subplot[0].set_xlim(0,22.32)
subplot[0].set_ylim(1.25,1.65)
subplot[1].set_xlim(22.34,44.66)
subplot[1].set_ylim(1.25,1.65)
subplot[2].set_xlim(44.66,66.98)
subplot[2].set_ylim(1.25,1.65)
subplot[1].set_ylim(1.25,1.65)
subplot[2].set_ylim(1.25,1.65)
subplot[0].legend(fontsize = 13, loc = 'upper left')






print('\n')
print('The new NEW method of computing # of events, stolen from the LCA')

def Hilton_LCA_AD_Leo(WTTS_M_flux, quiesence, sigma_1, time,announce):
    three_sigma_list_temp = []
    five_sigma_list_temp = []
    observed_data = WTTS_M_flux
    sigma_list = (observed_data-quiesence)/sigma_1
    for x in range(len(sigma_list)):
        sigma_list[x] = round(sigma_list[x],4)
    potential_flare_indecies = np.where(sigma_list > 2.5)[0]
    # potential_flare_indecies = [1,2,3,5,6,7,9]
    # potential_flare_indecies = [1,2,3,5,6,7,9,10,15,19,23,24]
    flares_list = []
    split_index = []
    for x in range(len(potential_flare_indecies)-1):
        if potential_flare_indecies[x+1]-potential_flare_indecies[x] > 1: split_index.append(x)
    if len(split_index) > 0:
        num_of_individ_flares = 'list_of_lists'
        for x in range(len(split_index)+1):
            # print(x)
            if x == 0:
                flares_list.append(np.array(potential_flare_indecies[0:split_index[0]+1]))
            elif x == len(split_index):
                  flares_list.append(np.array(potential_flare_indecies[split_index[-1]+1:]))
            else:
                flares_list.append(np.array(potential_flare_indecies[split_index[x-1]+1:split_index[x]+1]))
    elif len(split_index) == 0 and len(potential_flare_indecies) > 0:
        num_of_individ_flares = 'list_of_one'
        flares_list.append(potential_flare_indecies)
    else: num_of_individ_flares = 'list_of_none'
    if num_of_individ_flares == 'list_of_lists':
        flares_list[:] = [sublist for sublist in flares_list if any(3 <= sigma_list[x] for x in sublist) and  len(sublist) >= 3]
        #appending and decaring 3 sigma flares
        if len(flares_list) > 0:
            for x in flares_list:
                # print(x)
                loc_of_flare_peak = np.where(sigma_list == np.max(sigma_list[x]))[0]
                time_of_flare_peak = time[loc_of_flare_peak]
                three_sigma_list_temp.append(time_of_flare_peak)
                if announce == 'on':
                    print('a3 sigma flare detected at %s' % time_of_flare_peak)
        #deleting flares which don't breach the 5sigma level
        flares_list[:] = [sublist for sublist in flares_list if any(5 <= sigma_list[x] for x in sublist)]            
        #appending and decaring 5 sigma flares
        if len(flares_list) > 0:                
            for x in flares_list:
                loc_of_flare_peak = np.where(sigma_list == np.max(sigma_list[x]))[0]
                time_of_flare_peak = time[loc_of_flare_peak]
                five_sigma_list_temp.append(time_of_flare_peak)
                if announce == 'on':
                    print('b5 sigma flare detected at %s' % time_of_flare_peak)
                    # print('b5 sigma flare detected in %s at %s !!!!! img # %s' % (folder_list[i],time_of_flare_peak, i))
                
    if num_of_individ_flares == 'list_of_one' and len(flares_list[0]) >= 3:
        max_sig = np.max(sigma_list[tuple(flares_list)])
        if max_sig > 3:                               
            loc_of_flare_peak = np.where(sigma_list == max_sig)[0]
            time_of_flare_peak = time[loc_of_flare_peak]
            three_sigma_list_temp.append(time_of_flare_peak)
            if announce == 'on':
                print('c3 sigma flare detected at %s' % time_of_flare_peak)
                # print('c3 sigma flare detected in %s at %s ! img # %s' % (folder_list[i],time_of_flare_peak, i))
        if max_sig > 5:
            loc_of_flare_peak = np.where(sigma_list == max_sig)[0]
            time_of_flare_peak = time[loc_of_flare_peak]
            five_sigma_list_temp.append(time_of_flare_peak)
            if announce == 'on':
                print('d5 sigma flare detected at %s' % time_of_flare_peak)
                # print('d5 sigma flare detected in %s at %s !!!!! img # %s' % (folder_list[i],time_of_flare_peak, i))



    return three_sigma_list_temp, five_sigma_list_temp




# #######testing block
# three_sigma_list_temp = []
# five_sigma_list_temp = []
Fq_wtts_m = np.mean(ADLEO_TWA7_si4_flux[1000:1050])  #corresponding time is 20 to 21 ks
observed_data = ADLEO_TWA7_si4_flux - np.sqrt(ADLEO_TWA7_si4_error**2+popt_si4_high[2]**2)
quiesence = Fq_wtts_m
sigma_list = (observed_data-quiesence)/popt_si4_high[2]
# for x in range(len(sigma_list)):
#     sigma_list[x] = round(sigma_list[x],4)
# potential_flare_indecies = np.where(sigma_list > 2.5)[0]
# # potential_flare_indecies = [1,2,3,5,6,7,9]
# # potential_flare_indecies = [1,2,3,5,6,7,9,10,15,19,23,24]
# flares_list = []
# split_index = []
# for x in range(len(potential_flare_indecies)-1):
#     if potential_flare_indecies[x+1]-potential_flare_indecies[x] > 1: split_index.append(x)
# if len(split_index) > 0:
#     num_of_individ_flares = 'list_of_lists'
#     for x in range(len(split_index)+1):
#         # print(x)
#         if x == 0:
#             flares_list.append(np.array(potential_flare_indecies[0:split_index[0]+1]))
#         elif x == len(split_index):
#               flares_list.append(np.array(potential_flare_indecies[split_index[-1]+1:]))
#         else:
#             flares_list.append(np.array(potential_flare_indecies[split_index[x-1]+1:split_index[x]+1]))
# elif len(split_index) == 0 and len(potential_flare_indecies) > 0:
#     num_of_individ_flares = 'list_of_one'
#     flares_list.append(potential_flare_indecies)
# else: num_of_individ_flares = 'list_of_none'
# if num_of_individ_flares == 'list_of_lists':
#     flares_list[:] = [sublist for sublist in flares_list if any(3 <= sigma_list[x] for x in sublist) and  len(sublist) >= 3]
#     #appending and decaring 3 sigma flares
#     if len(flares_list) > 0:
#         for x in flares_list:
#             # print(x)
#             loc_of_flare_peak = np.where(sigma_list == np.max(sigma_list[x]))[0]
#             time_of_flare_peak = master_time[loc_of_flare_peak]
#             three_sigma_list_temp.append(time_of_flare_peak)
#             # if announce == 'on':
#             print('a3 sigma flare detected at %s' % time_of_flare_peak)
#     #deleting flares which don't breach the 5sigma level
#     flares_list[:] = [sublist for sublist in flares_list if any(5 <= sigma_list[x] for x in sublist)]            
#     #appending and decaring 5 sigma flares
#     if len(flares_list) > 0:                
#         for x in flares_list:
#             loc_of_flare_peak = np.where(sigma_list == np.max(sigma_list[x]))[0]
#             time_of_flare_peak = master_time[loc_of_flare_peak]
#             five_sigma_list_temp.append(time_of_flare_peak)
#             # if announce == 'on':
#             print('b5 sigma flare detected at %s' % time_of_flare_peak)
#                 # print('b5 sigma flare detected in %s at %s !!!!! img # %s' % (folder_list[i],time_of_flare_peak, i))
            
# if num_of_individ_flares == 'list_of_one' and len(flares_list[0]) >= 3:
#     max_sig = np.max(sigma_list[tuple(flares_list)])
#     if max_sig > 3:                               
#         loc_of_flare_peak = np.where(sigma_list == max_sig)[0]
#         time_of_flare_peak = master_time[loc_of_flare_peak]
#         three_sigma_list_temp.append(time_of_flare_peak)
#         # if announce == 'on':
#         print('c3 sigma flare detected at %s' % time_of_flare_peak)
#             # print('c3 sigma flare detected in %s at %s ! img # %s' % (folder_list[i],time_of_flare_peak, i))
#     if max_sig > 5:
#         loc_of_flare_peak = np.where(sigma_list == max_sig)[0]
#         time_of_flare_peak = master_time[loc_of_flare_peak]
#         five_sigma_list_temp.append(time_of_flare_peak)
#         # if announce == 'on':
#         print('d5 sigma flare detected at %s' % time_of_flare_peak)
#             # print('d5 sigma flare detected in %s at %s !!!!! img # %s' % (folder_list[i],time_of_flare_peak, i))

# ####################


three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
species_name = 'Si IV'
print('\n')
print('For %s:' % species_name)
a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_si4_flux, np.mean(ADLEO_TWA7_si4_flux[1000:1050]), popt_si4_high[2], master_time, announce = 'off')
for x in a:
    three_sigma_list.append(x)
for x in b:
    five_sigma_list.append(x)
   
a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_si4_flux+np.sqrt(ADLEO_TWA7_si4_error**2+popt_si4_high[2]**2), np.mean(ADLEO_TWA7_si4_flux[1000:1050]), popt_si4_high[2], master_time, announce = 'on')
for x in a:
    three_sigma_list_err_up.append(x)
for x in b:
    five_sigma_list_err_up.append(x)            

a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_si4_flux-np.sqrt(ADLEO_TWA7_si4_error**2+popt_si4_high[2]**2), np.mean(ADLEO_TWA7_si4_flux[1000:1050]), popt_si4_high[2], master_time, announce = 'off')
for x in a:
    three_sigma_list_err_down.append(x)
for x in b:
    five_sigma_list_err_down.append(x)


temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))



three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
species_name = 'C IV'
print('\n')
print('For %s:' % species_name)

a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_c4_flux, np.mean(ADLEO_TWA7_c4_flux[1000:1050]), popt_c4_high[2], master_time, announce = 'off')
for x in a:
    three_sigma_list.append(x)
for x in b:
    five_sigma_list.append(x)
   
a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_c4_flux+np.sqrt(ADLEO_TWA7_c4_error**2+popt_c4_high[2]**2), np.mean(ADLEO_TWA7_c4_flux[1000:1050]), popt_c4_high[2], master_time, announce = 'on')
for x in a:
    three_sigma_list_err_up.append(x)
for x in b:
    five_sigma_list_err_up.append(x)            

a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_c4_flux-np.sqrt(ADLEO_TWA7_c4_error**2+popt_c4_high[2]**2), np.mean(ADLEO_TWA7_c4_flux[1000:1050]), popt_c4_high[2], master_time, announce = 'off')
for x in a:
    three_sigma_list_err_down.append(x)
for x in b:
    five_sigma_list_err_down.append(x)


temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))


three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
species_name = 'He II'
print('\n')
print('For %s:' % species_name)
a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_He2_flux, np.mean(ADLEO_TWA7_He2_flux[1000:1050]), popt_He2_high[2], master_time, announce = 'off')
for x in a:
    three_sigma_list.append(x)
for x in b:
    five_sigma_list.append(x)
   
a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_He2_flux+np.sqrt(ADLEO_TWA7_He2_error**2+popt_He2_high[2]**2), np.mean(ADLEO_TWA7_He2_flux[1000:1050]), popt_He2_high[2], master_time, announce = 'on')
for x in a:
    three_sigma_list_err_up.append(x)
for x in b:
    five_sigma_list_err_up.append(x)            

a,b = Hilton_LCA_AD_Leo(ADLEO_TWA7_He2_flux-np.sqrt(ADLEO_TWA7_He2_error**2+popt_He2_high[2]**2), np.mean(ADLEO_TWA7_He2_flux[1000:1050]), popt_He2_high[2], master_time, announce = 'off')
for x in a:
    three_sigma_list_err_down.append(x)
for x in b:
    five_sigma_list_err_down.append(x)


temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))