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
    
    # plt.figure()
    # plt.plot(wave,flux)
    
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
master_si4_flux = master_si4_flux*10**12    
master_si4_error = master_si4_error*10**12

master_c4_flux = master_c4_flux*10**12    
master_c4_error = master_c4_error*10**12 

master_He2_flux = master_He2_flux*10**12    
master_He2_error = master_He2_error*10**12 

master_time = np.asarray(master_time); change_in_obs_flag = np.asarray(change_in_obs_flag);
master_time = master_time/1000; change_in_obs_flag = change_in_obs_flag/1000


#whole light curve plot
fig,subplot = plt.subplots(3, 1)
filename_counter = 0

for i in change_in_obs_flag:
    if i < master_time[1117]:
        subplot[0].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[0].text(x=i,y=1.5,s=filename_list[filename_counter][0:9], color = 'green', alpha = 0.75)
        subplot[0].text(x=i,y=1.2,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75)
        filename_counter+=1
subplot[0].plot(master_time[0:1117], master_si4_flux[0:1117], color = 'black', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], master_si4_flux[0:1117]-master_si4_error[0:1117],master_si4_flux[0:1117]+ master_si4_error[0:1117])
subplot[0].legend(fontsize = 16, loc = 'upper left')


for i in change_in_obs_flag:
    if i > master_time[1117] and i < master_time[2234]:
        subplot[1].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[1].text(x=i,y=8,s=filename_list[filename_counter][0:9], color = 'green', alpha = 0.75)
        subplot[1].text(x=i,y=7,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75)
        filename_counter+=1
subplot[1].plot(master_time[1117:2234], master_si4_flux[1117:2234], color = 'black', linewidth = 1)
subplot[1].fill_between(master_time[1117:2234], master_si4_flux[1117:2234]-master_si4_error[1117:2234],master_si4_flux[1117:2234]+ master_si4_error[1117:2234])


for i in change_in_obs_flag:
    if i > master_time[2234] and i < master_time[3349]:
        subplot[2].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[2].text(x=i,y=0.6,s=filename_list[filename_counter][0:9], color = 'green', alpha = 0.75)
        subplot[2].text(x=i,y=0.53,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75)
        filename_counter+=1
subplot[2].plot(master_time[2234:], master_si4_flux[2234:], color = 'black', linewidth = 1)
subplot[2].fill_between(master_time[2234:], master_si4_flux[2234:]-master_si4_error[2234:],master_si4_flux[2234:]+ master_si4_error[2234:])

# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Flux  ($\dfrac{erg}{cm^2s}$) *$10^{-12}$', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, '18.6hrs of STIS E140M Observations of AD Leo ($\Delta T = 20s$)', ha='center', va='center', fontsize = 18)
 
    

     





subplot0_HST_ID_fixes_list = [4.9,10.3,15.2,30.9,39]#manually selecting the change_in_obs_flag 's
                                  #that correspond with HST ID placement that needs
                                  #to be fixed

    
    
#whole light curve plot, y axis 0:1
fig,subplot = plt.subplots(3, 1)
filename_counter = 0

for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i < master_time[1117] and flag == False:
        subplot[0].axvline(x=i, color = 'green', alpha = 0.25, linewidth = 2)
        subplot[0].text(x=i,y=0.65,s=filename_list[filename_counter][0:9], color = 'green', alpha = 1)
        subplot[0].text(x=i,y=0.58,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 1)
    if i < master_time[1117]: filename_counter+=1
subplot[0].plot(master_time[0:1117], master_si4_flux[0:1117], color = 'black', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], master_si4_flux[0:1117]-master_si4_error[0:1117],master_si4_flux[0:1117]+ master_si4_error[0:1117])
subplot[0].legend(fontsize = 16, loc = 'upper left')
subplot[0].set_ylim(0,1)
#manually adding f/q ratios to large flares
subplot[0].text(x=13.07,y=0.9,s=r'$\dfrac{Flux}{Quiescence}=17$', color = 'magenta', alpha = 1)
subplot[0].arrow(x=15,y=0.93,dx=0.375, dy=0.03, color = 'magenta', alpha = 1, head_width = .034)
subplot[0].arrow(x=15,y=0.89,dx=0.54, dy=0.03, color = 'magenta', alpha = 1, head_width = .034)
#manually fixing HST IDs that need special placement
subplot[0].axvline(x=4.9, color = 'green', alpha = 0.25, linewidth = 4)
subplot[0].text(x=4.9725,y=0.9,s=filename_list[2][0:9], color = 'green', alpha = 0.75)
subplot[0].text(x=4.9725,y=0.83,s='MJD: ' + str(round(MJD_list[2],2)), color = 'green', alpha = 0.75)
subplot[0].axvline(x=10.3, color = 'green', alpha = 0.25, linewidth = 2)
subplot[0].text(x=10.3,y=0.9,s=filename_list[4][0:9], color = 'green', alpha = 0.75)
subplot[0].text(x=10.3,y=0.83,s='MJD: ' + str(round(MJD_list[4],2)), color = 'green', alpha = 0.75)
subplot[0].axvline(x=15.2, color = 'green', alpha = 0.25, linewidth = 2)
subplot[0].text(x=15.747,y=0.9,s=filename_list[6][0:9], color = 'green', alpha = 0.75)
subplot[0].text(x=15.747,y=0.83,s='MJD: ' + str(round(MJD_list[6],2)), color = 'green', alpha = 0.75)


for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > master_time[1117] and i < master_time[2234] and flag == False:
        subplot[1].axvline(x=i, color = 'green', alpha = 0.25, linewidth = 2)
        subplot[1].text(x=i,y=0.9,s=filename_list[filename_counter][0:9], color = 'green', alpha = 1)
        subplot[1].text(x=i,y=0.83,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 1)
    if i > master_time[1117] and i < master_time[2234]: filename_counter+=1
subplot[1].plot(master_time[1117:2234], master_si4_flux[1117:2234], color = 'black', linewidth = 1)
subplot[1].fill_between(master_time[1117:2234], master_si4_flux[1117:2234]-master_si4_error[1117:2234],master_si4_flux[1117:2234]+ master_si4_error[1117:2234])
subplot[1].set_ylim(0,1)
#manually adding f/q ratios to large flares
subplot[1].text(x=24.15,y=0.7,s='31', color = 'magenta', alpha = 1)
subplot[1].text(x=30.07,y=0.9,s='15', color = 'magenta', alpha = 1)
subplot[1].text(x=31.3,y=0.9,s='8', color = 'magenta', alpha = 1)
subplot[1].text(x=35.2,y=0.9,s='9', color = 'magenta', alpha = 1)
subplot[1].text(x=39.41,y=0.7,s='69', color = 'magenta', alpha = 1)
#manually fixing HST IDs that need special placement
subplot[1].axvline(x=30.9, color = 'green', alpha = 0.25, linewidth = 2)
subplot[1].text(x=31.73,y=0.9,s=filename_list[12][0:9], color = 'green', alpha = 1)
subplot[1].text(x=31.73,y=0.83,s='MJD: ' + str(round(MJD_list[12],2)), color = 'green', alpha = 0.75)
subplot[1].axvline(x=39, color = 'green', alpha = 0.25, linewidth = 2)
subplot[1].text(x=39.02,y=0.9,s=filename_list[15][0:9], color = 'green', alpha = 1)
subplot[1].text(x=39.02,y=0.83,s='MJD: ' + str(round(MJD_list[15],2)), color = 'green', alpha = 0.75)


for i in change_in_obs_flag:
    if i > master_time[2234] and i < master_time[3349]:
        subplot[2].axvline(x=i, color = 'green', alpha = 0.25, linewidth = 2)
        subplot[2].text(x=i,y=0.9,s=filename_list[filename_counter][0:9], color = 'green', alpha = 1)
        subplot[2].text(x=i,y=0.83,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 1)
        filename_counter+=1
subplot[2].plot(master_time[2234:], master_si4_flux[2234:], color = 'black', linewidth = 1)
subplot[2].fill_between(master_time[2234:], master_si4_flux[2234:]-master_si4_error[2234:],master_si4_flux[2234:]+ master_si4_error[2234:])
subplot[2].set_ylim(0,1)
# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.095, 0.5, r'Flux  ($erg \ cm^{-2} s^{-1}*10^{-12}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'STIS E140M Observations of AD Leo ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 








###########################


#mixing all of the ion species in one plot
fig,subplot = plt.subplots(1, 1)
filename_counter = 0

for i in change_in_obs_flag:
        subplot.axvline(x=i, color = 'black', alpha = 0.25, linewidth = 2)
        # subplot.text(x=i,y=1.2e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        # subplot.text(x=i+.23,y=1.2e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        # filename_counter+=1

subplot.plot(master_time, master_si4_flux, color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot.fill_between(master_time, master_si4_flux-master_si4_error,master_si4_flux+ master_si4_error, color = 'green',alpha = 0.25)    
subplot.plot(master_time, master_c4_flux, color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot.fill_between(master_time, master_c4_flux-master_c4_error,master_c4_flux+ master_c4_error,color = 'blue',alpha = 0.25)
subplot.plot(master_time, master_He2_flux, color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot.fill_between(master_time, master_He2_flux-master_He2_error,master_He2_flux+ master_He2_error, color = 'red',alpha = 0.25)


leg = subplot.legend(fontsize = 16, loc = 'upper right', markerscale =20)
leg.get_lines()[0].set_linewidth(3)
leg.get_lines()[1].set_linewidth(3)
leg.get_lines()[2].set_linewidth(3)
# subplot[0].set_ylim(-1e-12,4.5e-12)

# labels = [item.get_text() for item in subplot[0].get_yticklabels()]
# labels[1] = 0;labels[2] = 1;labels[3] = 2;labels[4] = 3;labels[5] = 4;
# subplot[0].set_yticklabels(labels)


fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.09, 0.5, r'Flux ($erg*cm^{-2}*s^{-1}*10^{-12}$)$', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'COS G160 Observations of TWA-7 ($\Delta t = %ss$, 48.6min Total)' % time_resolution, ha='center', va='center', fontsize = 18)




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
subplot[0].plot(master_time[0:1117], master_si4_flux[0:1117], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1117], master_si4_flux[0:1117]-master_si4_error[0:1117],master_si4_flux[0:1117]+ master_si4_error[0:1117],color = 'green', alpha = 0.1)

subplot[0].plot(master_time[0:1117], master_c4_flux[0:1117], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1117], master_c4_flux[0:1117]-master_c4_error[0:1117],master_c4_flux[0:1117]+ master_c4_error[0:1117],color = 'blue', alpha = 0.1)

subplot[0].plot(master_time[0:1117], master_He2_flux[0:1117], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[0].fill_between(master_time[0:1117], master_He2_flux[0:1117]-master_He2_error[0:1117],master_He2_flux[0:1117]+ master_He2_error[0:1117],color = 'red', alpha = 0.1)

subplot[0].legend(fontsize = 16, loc = 'upper left')
subplot[0].set_ylim(0,1)





subplot[1].plot(master_time[1117:2234], master_si4_flux[1117:2234], color = 'black', linewidth = 1)
subplot[1].fill_between(master_time[1117:2234], master_si4_flux[1117:2234]-master_si4_error[1117:2234],master_si4_flux[1117:2234]+ master_si4_error[1117:2234])
subplot[1].set_ylim(0,1)
#manually adding f/q ratios to large flares
subplot[1].text(x=24.15,y=0.7,s='31', color = 'magenta', alpha = 1)
subplot[1].text(x=30.07,y=0.9,s='15', color = 'magenta', alpha = 1)
subplot[1].text(x=31.3,y=0.9,s='8', color = 'magenta', alpha = 1)
subplot[1].text(x=35.2,y=0.9,s='9', color = 'magenta', alpha = 1)
subplot[1].text(x=39.41,y=0.7,s='69', color = 'magenta', alpha = 1)
#manually fixing HST IDs that need special placement
subplot[1].axvline(x=30.9, color = 'green', alpha = 0.25, linewidth = 2)
subplot[1].text(x=31.73,y=0.9,s=filename_list[12][0:9], color = 'green', alpha = 0.75)
subplot[1].text(x=31.73,y=0.83,s='MJD: ' + str(round(MJD_list[12],2)), color = 'green', alpha = 0.75)
subplot[1].axvline(x=39, color = 'green', alpha = 0.25, linewidth = 2)
subplot[1].text(x=39.02,y=0.9,s=filename_list[15][0:9], color = 'green', alpha = 0.75)
subplot[1].text(x=39.02,y=0.83,s='MJD: ' + str(round(MJD_list[15],2)), color = 'green', alpha = 0.75)


for i in change_in_obs_flag:
    if i > master_time[2234] and i < master_time[3349]:
        subplot[2].axvline(x=i, color = 'green', alpha = 0.25, linewidth = 2)
        subplot[2].text(x=i,y=0.9,s=filename_list[filename_counter][0:9], color = 'green', alpha = 0.75)
        subplot[2].text(x=i,y=0.83,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75)
        filename_counter+=1
subplot[2].plot(master_time[2234:], master_si4_flux[2234:], color = 'black', linewidth = 1)
subplot[2].fill_between(master_time[2234:], master_si4_flux[2234:]-master_si4_error[2234:],master_si4_flux[2234:]+ master_si4_error[2234:])
subplot[2].set_ylim(0,1)
# Set common labels
fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Flux  ($\dfrac{erg}{cm^2s}$) *$10^{-12}$', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'STIS E140M Observations of AD Leo ($\Delta t = %ss$, 18.6hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)
 

