# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 12:12:35 2020

@author: pahi9557

NOTE! Any changes to line locations must be done in two places (total + subdivided x1d's')

Time resolution is hard-coded
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

# def onclick(event):
#     global ix, iy, coords_x, coords_y, terminator
#     ix, iy = event.xdata, event.ydata
#     coords_x.append(ix)
#     coords_y.append(iy)
#     print('time:  '+ str(ix))
#     if len(coords_x) == 2:
#         fig.canvas.mpl_disconnect(fig)
#         plt.close()
#         coords_x = np.sort(coords_x)
#         terminator = 1
#     return

def subtract_quies(species_count_rate, quiescent_time_index):
    mean_quies = np.mean(species_count_rate[quiescent_time_index])
    subtracted_counts = species_count_rate - mean_quies
    return subtracted_counts

#function to find all tag files in directory
def get_file_names_with_strings(str_list):
    full_list = os.listdir(r"C:\Windows\System32\astroconda\astroconda_TW_Hya\COS_data\out20")
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]

    return final_list

def get_file_names_with_strings_PH(str_list,directory_path):
    # full_list = os.listdir(r'C:\Windows\System32\astroconda\AD_Leo_E140M')
    full_list = os.listdir(directory_path)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]

    return final_list

#defining the time resolution of this file
time_resolution = 20
#getting the name of all folders in /out, each one represents a full observation,
#the files inside are subframes
folder_list = get_file_names_with_strings(['l']) #all of the files
# This block of code was used to verify that the files I read in are in MJD
# sorted order
MJD_list = [] #defining empty MJD list, to be filled in 'for loop'
grating_list = []
for i in folder_list:
    filename = i + '-total-x1d.fits'
    # print(filename)
    hdul = fits.open(r"C:\Windows\System32\astroconda\astroconda_TW_Hya\COS_data\out20\%s\%s" % (i,filename))
    TEXPSTRT = hdul[1].header['EXPSTART'] #start time (MJD) of 1st exposure in file
    grating_list.append(hdul[0].header['OPT_ELEM'])
    MJD_list.append(TEXPSTRT)
    hdul.close()
sort_key = np.argsort(MJD_list) # this is the order that filename_list should be in
temp = []   # reparing empty variable to be made into the new filename_list                             
for i in sort_key: temp.append(folder_list[i])
folder_list = temp                   
MJD_list = np.asarray(MJD_list); MJD_list = MJD_list[sort_key]
del temp; del sort_key; # getting rid of finsihed variables   


temp = folder_list.index('le9d1cdiq'); del folder_list[temp]
temp = folder_list.index('le9d1cdkq'); del folder_list[temp]
temp = folder_list.index('le9d1cdmq'); del folder_list[temp]
temp = folder_list.index('le9d1cdoq'); del folder_list[temp]

temp = folder_list.index('le9d1kdfq'); del folder_list[temp]
temp = folder_list.index('le9d1kdhq'); del folder_list[temp]
temp = folder_list.index('le9d1kdjq'); del folder_list[temp]
temp = folder_list.index('le9d1kdmq'); del folder_list[temp]
del temp;


#preparing master lists
master_continuum_flux=np.asarray([]); master_continuum_counts=np.asarray([]); master_continuum_error=np.asarray([]);
master_si4_flux=np.asarray([]); master_si4_counts=np.asarray([]); master_si4_error=np.asarray([]);
master_c2_counts=np.asarray([]); master_c2_error=np.asarray([])
master_c3_counts=np.asarray([]); master_c3_error=np.asarray([])
master_c4_flux=np.asarray([]); master_c4_counts=np.asarray([]); master_c4_error=np.asarray([])
master_He2_flux=np.asarray([]); master_He2_counts=np.asarray([]); master_He2_error=np.asarray([]);
master_contin3_flux=np.asarray([]); master_contin3_counts=np.asarray([]); master_contin3_error=np.asarray([]);
master_contin4_flux=np.asarray([]); master_contin4_counts=np.asarray([]); master_contin4_error=np.asarray([]);
master_contin5_flux=np.asarray([]); master_contin5_counts=np.asarray([]); master_contin5_error=np.asarray([]);
master_H2_flux=np.asarray([]); master_H2_counts=np.asarray([]); master_H2_error=np.asarray([]);
master_commonBP_flux=np.asarray([]); master_commonBP_counts=np.asarray([]); master_commonBP_error=np.asarray([]);


master_time=[]; time = 0
change_in_obs_flag = [] 
cenwave_list = []

folder_count = 1
for i in tqdm(folder_list): # this loop iterates through each folder, a nested loop
#iterates through subframes and generates lightcurve
    print(' ');print(' ');
    print('Now working on %s of %s G160 folders' % (folder_count, len(folder_list)))
    change_in_obs_flag.append(time)
    folder = i
    directory_path = r"C:\Windows\System32\astroconda\astroconda_TW_Hya\COS_data\out20\%s" % folder
    filelist = get_file_names_with_strings_PH(['split'], directory_path)
    num_of_obs = int(len(filelist)/2)
    for j in range(num_of_obs): # this loop iterates through each sub-divided image, but starts with
    # the total file via an if statement
        #generates lightcurve
        split_counter = j + 1

        #need to do the total file only once for each folder, i.e. on split_1
        if split_counter == 1:
            #locating the corresponding 'total' file for flux callibration
            filename_tot = i + "-total-x1d.fits"
            hdul_tot  = fits.open(r"C:\Windows\System32\astroconda\astroconda_TW_Hya\COS_data\out20\%s\%s" % (i, filename_tot))

            data = hdul_tot[1].data
            cenwave_list.append(hdul_tot[0].header['CENWAVE'])
            FUVa = data[0]; FUVb = data[1];
            wave_a = FUVa['wavelength']; wave_b = FUVb['wavelength'];
            flux_a = FUVa['flux']; flux_b = FUVb['flux'];
            count_a = FUVa['net']; count_b = FUVb['net'];
            

            tot_x1d_wav_array = np.concatenate((wave_b, wave_a))
            tot_x1d_flux_array = np.concatenate((flux_b, flux_a))
            tot_x1d_count_array = np.concatenate((count_b, count_a))
            
            wave = tot_x1d_wav_array; wave = wave.flatten()
            flux = tot_x1d_flux_array; flux = flux.flatten()
            counts = tot_x1d_count_array; counts = counts.flatten()
            
            sort_key = np.argsort(wave)
            wave = wave[sort_key]; flux = flux[sort_key]; counts = counts[sort_key]
            
            # plt.figure()
            # plt.plot(wave,flux)
            
            
             ################ indentifying line locations############
            c3_line_location = np.where((wave > 1174) & (wave < 1176.5))
            n5_line_location = np.where((wave > 1238) & (wave < 1243.3))
            c4_line_location = np.where((wave > 1547.7) & (wave < 1553.35))
            He2_line_location = np.where((wave > 1639.6) & (wave < 1641.5))
            
            c2_line_location = np.where((wave > 1334) & (wave < 1336.5))
            si2_line_location = np.where((wave > 1264.34) & (wave < 1265.5))
            si3_line_location = np.where((wave > 1205.9) & (wave < 1207.3))
            o1_line_location = np.where((wave > 1304.6) & (wave < 1306.5))
            
            continuum_line_location = np.where((wave > 1339) & (wave < 1351))
            
            bandpass1_1 =  np.asarray(np.where((wave > 1439 ) & (wave < 1561))) #common b segment
            bandpass1_2 =  np.asarray(np.where((wave > 1635 ) & (wave < 1700))) #common a segment
            commonBP_line_location = []
            for j in range(np.size(bandpass1_1)):
                commonBP_line_location.append(bandpass1_1[0][j])
            for j in range(np.size(bandpass1_2)):
                commonBP_line_location.append(bandpass1_2[0][j])
            
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
            
            #now for continuum # 2 and the H2 lines in that region
            
            contin4_line_location = np.where((wave > 1644) & (wave < 1655))
            contin5_line_location = np.where((wave > 1680) & (wave < 1750))
            
            contin3_line_location =  []
            bandpass1 =  np.asarray(np.where((wave > 1436 ) & (wave < 1442)))
            bandpass2 =  np.asarray(np.where((wave > 1447 ) & (wave < 1452)))
            bandpass3 =  np.asarray(np.where((wave > 1468 ) & (wave < 1487)))
            bandpass4 =  np.asarray(np.where((wave > 1508 ) & (wave < 1513)))
            for j in range(np.size(bandpass1)):
                contin3_line_location.append(bandpass1[0][j])
            for j in range(np.size(bandpass2)):
                contin3_line_location.append(bandpass2[0][j])
            for j in range(np.size(bandpass3)):
                contin3_line_location.append(bandpass3[0][j])
            for j in range(np.size(bandpass4)):
                contin3_line_location.append(bandpass4[0][j])

            H2_line_location = []
            bandpass1 =  np.asarray(np.where((wave > 1445 ) & (wave < 1447)))
            bandpass2 =  np.asarray(np.where((wave > 1452 ) & (wave < 1455)))
            bandpass3 =  np.asarray(np.where((wave > 1462 ) & (wave < 1465)))
            bandpass4 =  np.asarray(np.where((wave > 1488 ) & (wave < 1490)))
            bandpass5 =  np.asarray(np.where((wave > 1499 ) & (wave < 1501)))
            bandpass6 =  np.asarray(np.where((wave > 1504 ) & (wave < 1506)))

            for j in range(np.size(bandpass1)):
                H2_line_location.append(bandpass1[0][j])
            for j in range(np.size(bandpass2)):
                H2_line_location.append(bandpass2[0][j])
            for j in range(np.size(bandpass3)):
                H2_line_location.append(bandpass3[0][j])
            for j in range(np.size(bandpass4)):
                H2_line_location.append(bandpass4[0][j])
            for j in range(np.size(bandpass5)):
                H2_line_location.append(bandpass5[0][j])
            for j in range(np.size(bandpass6)):
                H2_line_location.append(bandpass6[0][j])
            
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
            
            avg_contin3_line_flux = np.trapz(flux[contin3_line_location],wave[contin3_line_location])
            avg_contin4_line_flux = np.trapz(flux[contin4_line_location],wave[contin4_line_location])
            avg_contin5_line_flux = np.trapz(flux[contin5_line_location],wave[contin5_line_location])
            avg_H2_line_flux = np.trapz(flux[H2_line_location],wave[H2_line_location])
            avg_commonBP_line_flux = np.trapz(flux[commonBP_line_location],wave[commonBP_line_location])
            
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
            
            avg_contin3_line_counts = np.sum(counts[contin3_line_location])
            avg_contin4_line_counts = np.sum(counts[contin4_line_location])
            avg_contin5_line_counts = np.sum(counts[contin5_line_location])
            avg_H2_line_counts = np.sum(counts[H2_line_location])
            avg_commonBP_line_counts = np.sum(counts[commonBP_line_location])
            
            hdul_tot.close()
            

        
        filename_a = i + "_split_a_%s_x1d.fits" % split_counter
        filename_b = i + "_split_b_%s_x1d.fits" % split_counter
        hdula = fits.open(r"C:\Windows\System32\astroconda\astroconda_TW_Hya\COS_data\out20\%s\%s" % (i, filename_a))
        hdulb = fits.open(r"C:\Windows\System32\astroconda\astroconda_TW_Hya\COS_data\out20\%s\%s" % (i, filename_b))
        HST_ID = i #the first 9 characters of an HST filename is the ID
        TEXPTIME = hdula[1].header['EXPTIMEa'] #total exposure time in seconds
        TEXPTIME = int(TEXPTIME) #rounding to the nearest second
        # print('Total Exposure time =  %ss' % TEXPTIME)
        
        data_a = hdula[1].data; data_b = hdulb[1].data;
        FUVa = data_a[0]; FUVb = data_b[0];
        wave_a = FUVa['wavelength']; wave_b = FUVb['wavelength']; 
        flux_a = FUVa['flux']; flux_b = FUVb['flux']; 
        count_a = FUVa['net']; count_b = FUVb['net'];
        
        wave = np.concatenate((wave_b, wave_a))
        counts = np.concatenate((count_b, count_a))
        sort_key = np.argsort(wave)
        wave = wave[sort_key]; counts = counts[sort_key]

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
        contin3_flux=[]; contin3_flux_error=[]; contin3_counts=[]; contin3_error=[];
        contin3_quiescent_flux=[]; contin3_quiescent_flux_error=[]; contin3_quiescent_counts=[]; contin3_quiescent_error=[];
        contin4_flux=[]; contin4_flux_error=[]; contin4_counts=[]; contin4_error=[];
        contin5_flux=[]; contin5_flux_error=[]; contin5_counts=[]; contin5_error=[];
        H2_flux=[]; H2_flux_error=[]; H2_counts=[]; H2_error=[];
        commonBP_flux=[]; commonBP_flux_error=[]; commonBP_counts=[]; commonBP_error=[];
        bandpass1_flux=[]; bandpass1_flux_error=[]; bandpass1_counts=[]; bandpass1_error=[];
        bandpass2_flux=[]; bandpass2_flux_error=[]; bandpass2_counts=[]; bandpass2_error=[];
        
         ################ indentifying line locations############
        c3_line_location = np.where((wave > 1174) & (wave < 1176.5))
        n5_line_location = np.where((wave > 1238) & (wave < 1243.3))
        c4_line_location = np.where((wave > 1547.7) & (wave < 1553.35))
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
        
        
        #now for continuum # 2 and the H2 lines in that region
        
        contin4_line_location = np.where((wave > 1644) & (wave < 1655))
        contin5_line_location = np.where((wave > 1680) & (wave < 1750))
        
        bandpass1_1 =  np.asarray(np.where((wave > 1439 ) & (wave < 1561))) #common b segment
        bandpass1_2 =  np.asarray(np.where((wave > 1635 ) & (wave < 1700))) #common a segment
        commonBP_line_location = []
        for j in range(np.size(bandpass1_1)):
            commonBP_line_location.append(bandpass1_1[0][j])
        for j in range(np.size(bandpass1_2)):
            commonBP_line_location.append(bandpass1_2[0][j])
    
        contin3_line_location =  []
        bandpass1 =  np.asarray(np.where((wave > 1436 ) & (wave < 1442)))
        bandpass2 =  np.asarray(np.where((wave > 1447 ) & (wave < 1452)))
        bandpass3 =  np.asarray(np.where((wave > 1468 ) & (wave < 1487)))
        bandpass4 =  np.asarray(np.where((wave > 1508 ) & (wave < 1513)))
        for j in range(np.size(bandpass1)):
            contin3_line_location.append(bandpass1[0][j])
        for j in range(np.size(bandpass2)):
            contin3_line_location.append(bandpass2[0][j])
        for j in range(np.size(bandpass3)):
            contin3_line_location.append(bandpass3[0][j])
        for j in range(np.size(bandpass4)):
            contin3_line_location.append(bandpass4[0][j])

        H2_line_location = []
        bandpass1 =  np.asarray(np.where((wave > 1445 ) & (wave < 1447)))
        bandpass2 =  np.asarray(np.where((wave > 1452 ) & (wave < 1455)))
        bandpass3 =  np.asarray(np.where((wave > 1462 ) & (wave < 1465)))
        bandpass4 =  np.asarray(np.where((wave > 1488 ) & (wave < 1490)))
        bandpass5 =  np.asarray(np.where((wave > 1499 ) & (wave < 1501)))
        bandpass6 =  np.asarray(np.where((wave > 1504 ) & (wave < 1506)))

        for j in range(np.size(bandpass1)):
            H2_line_location.append(bandpass1[0][j])
        for j in range(np.size(bandpass2)):
            H2_line_location.append(bandpass2[0][j])
        for j in range(np.size(bandpass3)):
            H2_line_location.append(bandpass3[0][j])
        for j in range(np.size(bandpass4)):
            H2_line_location.append(bandpass4[0][j])
        for j in range(np.size(bandpass5)):
            H2_line_location.append(bandpass5[0][j])
        for j in range(np.size(bandpass6)):
            H2_line_location.append(bandpass6[0][j])
        
        ################################################################
        
        # computing count rates
        c2temp_counts = np.sum(counts[c2_line_location])
        c3temp_counts = np.sum(counts[c3_line_location])
        c4temp_counts = np.sum(counts[c4_line_location])
        
        continuum_temp_counts = np.sum(counts[continuum_line_location])
        si4_temp_counts = np.sum(counts[si4_line_location])
        He2_temp_counts = np.sum(counts[He2_line_location])
        
        contin3_temp_counts = np.sum(counts[contin3_line_location])
        contin4_temp_counts = np.sum(counts[contin4_line_location])
        contin5_temp_counts = np.sum(counts[contin5_line_location])
        
        #contin3 quiescent rate, take the mean instead of the sum
        contin3_temp_quiescent_counts = np.mean(counts[contin3_line_location])
        
        H2_temp_counts = np.sum(counts[H2_line_location]-contin3_temp_quiescent_counts)
        commonBP_temp_counts = np.sum(counts[commonBP_line_location])
        
        
        #appending count rates
        c2_counts.append(c2temp_counts); c2_error.append(np.sqrt(c2temp_counts*time_resolution)/time_resolution)
        c3_counts.append(c3temp_counts); c3_error.append(np.sqrt(c3temp_counts*time_resolution)/time_resolution)
        c4_counts.append(c4temp_counts); c4_error.append(np.sqrt(c4temp_counts*time_resolution)/time_resolution)
        
        continuum_counts.append(continuum_temp_counts); continuum_error.append(np.sqrt(continuum_temp_counts*time_resolution)/time_resolution)
        si4_counts.append(si4_temp_counts); si4_error.append(np.sqrt(si4_temp_counts*time_resolution)/time_resolution)
        He2_counts.append(He2_temp_counts); He2_error.append(np.sqrt(He2_temp_counts*time_resolution)/time_resolution)
        
        contin3_counts.append(contin3_temp_counts); contin3_error.append(np.sqrt(contin3_temp_counts*time_resolution)/time_resolution)
        contin4_counts.append(contin4_temp_counts); contin4_error.append(np.sqrt(contin4_temp_counts*time_resolution)/time_resolution)
        contin5_counts.append(contin5_temp_counts); contin5_error.append(np.sqrt(contin5_temp_counts*time_resolution)/time_resolution)
        H2_counts.append(H2_temp_counts); H2_error.append(np.sqrt(H2_temp_counts*time_resolution)/time_resolution)
        commonBP_counts.append(commonBP_temp_counts); commonBP_error.append(np.sqrt(commonBP_temp_counts*time_resolution)/time_resolution)
        
        
        #converting count lists to arrays
        c2_counts=np.asarray(c2_counts); c2_error=np.asarray(c2_error);
        c3_counts=np.asarray(c3_counts); c3_error=np.asarray(c3_error);
        c4_counts=np.asarray(c4_counts); c4_error=np.asarray(c4_error);
        
        continuum_counts=np.asarray(continuum_counts); continuum_error=np.asarray(continuum_error);
        si4_counts=np.asarray(si4_counts); si4_error=np.asarray(si4_error);
        He2_counts=np.asarray(He2_counts); He2_error=np.asarray(He2_error);
        
        contin3_counts=np.asarray(contin3_counts); contin3_error=np.asarray(contin3_error);
        contin4_counts=np.asarray(contin4_counts); contin4_error=np.asarray(contin4_error);
        contin5_counts=np.asarray(contin5_counts); contin5_error=np.asarray(contin5_error);
        H2_counts=np.asarray(H2_counts); H2_error=np.asarray(H2_error);
        commonBP_counts=np.asarray(commonBP_counts); commonBP_error=np.asarray(commonBP_error);
        
        
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
        
        fluxed_contin3 = contin3_counts*(avg_contin3_line_flux/avg_contin3_line_counts)
        fluxed_contin3_err = contin3_error*(avg_contin3_line_flux/avg_contin3_line_counts)
        fluxed_contin4 = contin4_counts*(avg_contin4_line_flux/avg_contin4_line_counts)
        fluxed_contin4_err = contin4_error*(avg_contin4_line_flux/avg_contin4_line_counts)
        fluxed_contin5 = contin5_counts*(avg_contin5_line_flux/avg_contin5_line_counts)
        fluxed_contin5_err = contin5_error*(avg_contin5_line_flux/avg_contin5_line_counts)
        fluxed_H2 = H2_counts*(avg_H2_line_flux/avg_H2_line_counts)
        fluxed_H2_err = H2_error*(avg_H2_line_flux/avg_H2_line_counts)
        fluxed_commonBP = commonBP_counts*(avg_commonBP_line_flux/avg_commonBP_line_counts)
        fluxed_commonBP_err = commonBP_error*(avg_commonBP_line_flux/avg_commonBP_line_counts)
        
    
        master_continuum_counts= np.append(master_continuum_counts,continuum_counts);
        master_continuum_flux=   np.append(master_continuum_flux,fluxed_continuum); 
        master_continuum_error=  np.append(master_continuum_error,fluxed_continuum_err); 
        
        master_c4_counts= np.append(master_c4_counts,c4_counts);
        master_c4_flux=   np.append(master_c4_flux,fluxed_c4); 
        master_c4_error=  np.append(master_c4_error,fluxed_c4_err); 
        
        master_si4_counts= np.append(master_si4_counts,si4_counts);
        master_si4_flux=   np.append(master_si4_flux,fluxed_si4); 
        master_si4_error=  np.append(master_si4_error,fluxed_si4_err); 
        
        master_He2_counts= np.append(master_He2_counts,He2_counts);
        master_He2_flux=   np.append(master_He2_flux,fluxed_He2); 
        master_He2_error=  np.append(master_He2_error,fluxed_He2_err); 
        
        master_contin3_counts= np.append(master_contin3_counts,contin3_counts);
        master_contin3_flux=   np.append(master_contin3_flux,fluxed_contin3); 
        master_contin3_error=  np.append(master_contin3_error,fluxed_contin3_err); 
        
        master_contin4_counts= np.append(master_contin4_counts,contin4_counts);
        master_contin4_flux=   np.append(master_contin4_flux,fluxed_contin4); 
        master_contin4_error=  np.append(master_contin4_error,fluxed_contin4_err); 
        
        master_contin5_counts= np.append(master_contin5_counts,contin5_counts);
        master_contin5_flux=   np.append(master_contin5_flux,fluxed_contin5); 
        master_contin5_error=  np.append(master_contin5_error,fluxed_contin5_err); 
        
        master_H2_counts= np.append(master_H2_counts,H2_counts);
        master_H2_flux=   np.append(master_H2_flux,fluxed_H2); 
        master_H2_error=  np.append(master_H2_error,fluxed_H2_err); 
    
        master_commonBP_counts= np.append(master_commonBP_counts,commonBP_counts);
        master_commonBP_flux=   np.append(master_commonBP_flux,fluxed_commonBP); 
        master_commonBP_error=  np.append(master_commonBP_error,fluxed_commonBP_err); 
        
         
        
        
        
        master_time.append(time)
        time = time + time_resolution
    folder_count += 1 
        
        
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
del c3_line_location;del c3temp_counts; del contin3_counts;del contin3_error;del contin3_flux;del contin3_flux_error;del contin3_line_location;
del avg_contin3_line_counts;del avg_contin3_line_flux; del contin4_counts;del contin4_error;del contin4_flux;del contin4_flux_error;del contin4_line_location;
del avg_contin4_line_counts;del avg_contin4_line_flux; del contin5_counts;del contin5_error;del contin5_flux;del contin5_flux_error;del contin5_line_location;
del avg_contin5_line_counts;del avg_contin5_line_flux; del H2_counts;del H2_error;del H2_flux;del H2_flux_error;del H2_line_location;
del avg_H2_line_counts;del avg_H2_line_flux;

 

#scaling
master_time = np.asarray(master_time); change_in_obs_flag = np.asarray(change_in_obs_flag);
master_time = master_time/1000; change_in_obs_flag = change_in_obs_flag/1000


subplot0_HST_ID_fixes_list = [5.28, 5.8, 6.32, 9.96, 10.48, 11,
                              16.2, 16.88, 17.56, 18.24, 18.92, 19.6,
                              20.28, 20.96, 21.64, 40.68]#manually selecting the change_in_obs_flag 's
                                  #that correspond with HST ID placement that needs
                                  #to be fixed


#whole si4 light curve plot
fig,subplot = plt.subplots(2, 1)
filename_counter = 0

for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i < master_time[1107] and flag == False:
        subplot[0].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[0].text(x=i,y=0.4e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        subplot[0].text(x=i+.23,y=0.4e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i < master_time[1107] and flag == True:
        subplot[0].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[0].text(x=i,y=1e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        subplot[0].text(x=i+.23,y=1e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    
subplot[0].plot(master_time[0:1107], master_si4_flux[0:1107], color = 'black', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1107], master_si4_flux[0:1107]-master_si4_error[0:1107],master_si4_flux[0:1107]+ master_si4_error[0:1107])
subplot[0].legend(fontsize = 16, loc = 'upper left')

for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > master_time[1106] and flag == False:
        subplot[1].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[1].text(x=i,y=1.36e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        subplot[1].text(x=i+.23,y=1.36e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i > master_time[1106] and flag == True:
        subplot[1].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[1].text(x=i,y=.6e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        subplot[1].text(x=i+.23,y=.6e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
subplot[1].plot(master_time[1107:], master_si4_flux[1107:], color = 'black', linewidth = 1)
subplot[1].fill_between(master_time[1107:], master_si4_flux[1107:]-master_si4_error[1107:],master_si4_flux[1107:]+ master_si4_error[1107:])

fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Flux  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'COS G160 Observations of TW Hydrae ($\Delta t = %ss$, 12.4hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)



#re-doing for c4 / He2
subplot0_HST_ID_fixes_list = [5.28, 5.8, 6.32, 9.96, 10.48, 11,
                              16.2, 16.88, 17.56, 18.24, 18.92, 19.6,
                              20.28, 20.96, 21.64, 40.68, 44.08, 43.4, 42.72, 42.04, 41.36,
                              40.0, 39.32, 38.64, 37.96, 37.28]#manually selecting the change_in_obs_flag 's
                                  #that correspond with HST ID placement that needs
                                  #to be fixed


#whole c4 light curve plot
fig,subplot = plt.subplots(2, 1)
filename_counter = 0

for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i < master_time[1107] and flag == False:
        subplot[0].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[0].text(x=i,y=1.2e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        subplot[0].text(x=i+.23,y=1.2e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i < master_time[1107] and flag == True:
        subplot[0].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[0].text(x=i,y=3e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        subplot[0].text(x=i+.23,y=3e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    
subplot[0].plot(master_time[0:1107], master_c4_flux[0:1107], color = 'black', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1107], master_c4_flux[0:1107]-master_c4_error[0:1107],master_c4_flux[0:1107]+ master_c4_error[0:1107])
subplot[0].legend(fontsize = 16, loc = 'upper left')

for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > master_time[1106] and flag == False:
        subplot[1].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[1].text(x=i,y=3e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        subplot[1].text(x=i+.23,y=3e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i > master_time[1106] and flag == True:
        subplot[1].axvline(x=i, color = 'green', alpha = 0.25)
        subplot[1].text(x=i,y=1.1e-12,s=folder_list[filename_counter][0:9], color = 'green', alpha = 0.75, rotation = 'vertical')
        subplot[1].text(x=i+.23,y=1.1e-12,s='MJD: ' + str(round(MJD_list[filename_counter],2)), color = 'green', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
subplot[1].plot(master_time[1107:], master_c4_flux[1107:], color = 'black', linewidth = 1)
subplot[1].fill_between(master_time[1107:], master_c4_flux[1107:]-master_c4_error[1107:],master_c4_flux[1107:]+ master_c4_error[1107:])
subplot[1].set_ylim(1e-12,4.5e-12)

fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.075, 0.5, r'Flux  ($\dfrac{erg}{cm^2s}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'COS G160 Observations of TW Hydrae ($\Delta t = %ss$, 12.4hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)







#mixing all of the ion species in one plot
fig,subplot = plt.subplots(3, 1, figsize = [19.3,9.5])
filename_counter = 0

for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i < master_time[1116] and flag == False:
        subplot[0].axvline(x=i, color = 'black', alpha = 0.5)
        subplot[0].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        subplot[0].text(x=i+.23,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i < master_time[1116] and flag == True:
        subplot[0].axvline(x=i, color = 'black', alpha = 0.5)
        subplot[0].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        subplot[0].text(x=i+.23,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1

subplot[0].plot(master_time[0:1116], 1e12*master_si4_flux[0:1116], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
subplot[0].fill_between(master_time[0:1116], 1e12*master_si4_flux[0:1116]-1e12*master_si4_error[0:1116],1e12*master_si4_flux[0:1116]+ 1e12*master_si4_error[0:1116], color = 'green',alpha = 0.25)    
subplot[0].plot(master_time[0:1116], 1e12*master_c4_flux[0:1116], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
subplot[0].fill_between(master_time[0:1116], 1e12*master_c4_flux[0:1116]-1e12*master_c4_error[0:1116],1e12*master_c4_flux[0:1116]+ 1e12*master_c4_error[0:1116],color = 'blue',alpha = 0.25)
subplot[0].plot(master_time[0:1116], 1e12*master_He2_flux[0:1116], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
subplot[0].fill_between(master_time[0:1116], 1e12*master_He2_flux[0:1116]-1e12*master_He2_error[0:1116],1e12*master_He2_flux[0:1116]+ 1e12*master_He2_error[0:1116], color = 'red',alpha = 0.25)
leg = subplot[0].legend(fontsize = 12, loc = 'upper right', markerscale =12)
leg.get_lines()[0].set_linewidth(6)
leg.get_lines()[1].set_linewidth(6)
leg.get_lines()[2].set_linewidth(6)
subplot[0].set_ylim(1e12*(-2.5e-12),1e12*5e-12)
subplot[0].set_xlim(-.25,22.5)

# labels = [item.get_text() for item in subplot[0].get_yticklabels()]
# labels[0] = ''; labels[1] = '';
# # labels[1] = 0;labels[2] = 1;labels[3] = 2;labels[4] = 3;labels[5] = 4;
# subplot[0].set_yticklabels(labels)


for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > master_time[1115] and i < master_time[2238] and flag == False:
        subplot[1].axvline(x=i, color = 'black', alpha = 0.5)
        subplot[1].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        subplot[1].text(x=i+.23,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i > master_time[1115] and i < master_time[2238] and flag == True:
        subplot[1].axvline(x=i, color = 'black', alpha = 0.5)
        subplot[1].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        subplot[1].text(x=i+.23,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1

subplot[1].plot(master_time[1116:2238], 1e12*master_si4_flux[1116:2238], color = 'green', linewidth = 1)
subplot[1].fill_between(master_time[1116:2238], 1e12*master_si4_flux[1116:2238]-1e12*master_si4_error[1116:2238],1e12*master_si4_flux[1116:2238]+ 1e12*master_si4_error[1116:2238], color = 'green',alpha = 0.25)
subplot[1].plot(master_time[1116:2238], 1e12*master_c4_flux[1116:2238], color = 'blue', linewidth = 1)
subplot[1].fill_between(master_time[1116:2238], 1e12*master_c4_flux[1116:2238]-1e12*master_c4_error[1116:2238],1e12*master_c4_flux[1116:2238]+ 1e12*master_c4_error[1116:2238],color = 'blue', alpha = 0.25)
subplot[1].plot(master_time[1116:2238], 1e12*master_He2_flux[1116:2238], color = 'red', linewidth = 1)
subplot[1].fill_between(master_time[1116:2238], 1e12*master_He2_flux[1116:2238]-1e12*master_He2_error[1116:2238],1e12*master_He2_flux[1116:2238]+ 1e12*master_He2_error[1116:2238], color = 'red',alpha = 0.25)
subplot[1].set_ylim(1e12*(-2.5e-12),1e12*5e-12)
subplot[1].set_xlim(22.1,45)

# labels = [item.get_text() for item in subplot[1].get_yticklabels()]
# labels[1] = 0;labels[2] = 1;labels[3] = 2;labels[4] = 3;labels[5] = 4;
# subplot[1].set_yticklabels(labels)


for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > master_time[2237] and flag == False:
        subplot[2].axvline(x=i, color = 'black', alpha = 0.5)
        subplot[2].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        subplot[2].text(x=i+.13,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i > master_time[2237] and flag == True:
        subplot[2].axvline(x=i, color = 'black', alpha = 0.5)
        subplot[2].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        subplot[2].text(x=i+.13,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1

subplot[2].plot(master_time[2238:], 1e12*master_si4_flux[2238:], color = 'green', linewidth = 1)
subplot[2].fill_between(master_time[2238:], 1e12*master_si4_flux[2238:]-1e12*master_si4_error[2238:],1e12*master_si4_flux[2238:]+ 1e12*master_si4_error[2238:], color = 'green',alpha = 0.25)
subplot[2].plot(master_time[2238:], 1e12*master_c4_flux[2238:], color = 'blue', linewidth = 1)
subplot[2].fill_between(master_time[2238:], 1e12*master_c4_flux[2238:]-1e12*master_c4_error[2238:],1e12*master_c4_flux[2238:]+ 1e12*master_c4_error[2238:],color = 'blue', alpha = 0.25)
subplot[2].plot(master_time[2238:], 1e12*master_He2_flux[2238:], color = 'red', linewidth = 1)
subplot[2].fill_between(master_time[2238:], 1e12*master_He2_flux[2238:]-1e12*master_He2_error[2238:],1e12*master_He2_flux[2238:]+ 1e12*master_He2_error[2238:], color = 'red',alpha = 0.25)
subplot[2].set_ylim(1e12*(-2.5e-12),1e12*5e-12)
subplot[2].set_xlim(44.65,56.85)


fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.095, 0.5, r'Flux  ($erg \ cm^{-2} s^{-1}*10^{-12}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'COS G160M Observations of TW Hydrae ($\Delta t = %ss$, 15.75hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)






temp_time = np.linspace(1, len(master_time[2238:]),len(master_time[2238:]))*20
filename_counter = 72
plt.figure()
# for i in change_in_obs_flag:
#     flag = i in subplot0_HST_ID_fixes_list
#     if i > master_time[2237] and flag == False:
#         plt.axvline(x=i, color = 'black', alpha = 0.5)
#         plt.text(x=i,y=-2.3e-12,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
#         plt.text(x=i+.13,y=-2.3e-12,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
#         filename_counter+=1
#     if i > master_time[2237] and flag == True:
#         plt.axvline(x=i, color = 'black', alpha = 0.5)
#         # plt.text(x=i,y=-2.3e-12,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
#         # plt.text(x=i+.13,y=-2.3e-12,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
#         # filename_counter+=1
temp_change_in_obs_flag = (change_in_obs_flag - 44.76)*1000
for i in temp_change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > -1:
        plt.axvline(x=i, color = 'black', alpha = 0.5)
        # plt.text(x=i,y=-2.3e-12,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        # plt.text(x=i+.13,y=-2.3e-12,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        # filename_counter+=1

# plt.plot(temp_time, master_si4_flux[2238:], color = 'green', linewidth = 1)
# plt.fill_between(temp_time, master_si4_flux[2238:]-master_si4_error[2238:],master_si4_flux[2238:]+ master_si4_error[2238:], color = 'green',alpha = 0.25)
plt.plot(temp_time, master_c4_flux[2238:]-3.28e-12, color = 'green', linewidth = 1, label = 'C IV')
plt.fill_between(temp_time, master_c4_flux[2238:]-master_c4_error[2238:]-3.28e-12,master_c4_flux[2238:]+ master_c4_error[2238:]-3.28e-12,color = 'green', alpha = 0.25)

plt.plot(temp_time, master_H2_flux[2238:]-1.23e-12, color = 'red', linewidth = 1, label = 'H2')
plt.fill_between(temp_time, master_H2_flux[2238:]-master_H2_error[2238:]-1.23e-12,master_H2_flux[2238:]+ master_H2_error[2238:]-1.23e-12, color = 'red',alpha = 0.25)

plt.plot(temp_time, master_contin3_flux[2238:]-2.06e-12, color = 'blue', linewidth = 1, label = 'Continuum in (1430 - 1520)$\AA$')
plt.fill_between(temp_time, master_contin3_flux[2238:]-master_contin3_error[2238:]-2.06e-12,master_contin3_flux[2238:]+ master_contin3_error[2238:]-2.06e-12, color = 'blue',alpha = 0.25)

plt.plot(temp_time, master_contin4_flux[2238:]+.43e-12, color = 'orange', linewidth = 1, label = 'Continuum in (1600 - 1655)$\AA$')
plt.fill_between(temp_time, master_contin4_flux[2238:]-master_contin4_error[2238:]+.43e-12,master_contin4_flux[2238:]+ master_contin4_error[2238:]+.43e-12, color = 'orange',alpha = 0.25)

plt.plot(temp_time, master_contin5_flux[2238:]-2.54e-12, color = 'magenta', linewidth = 1, label = 'Continuum in (1680 - 1750)$\AA$')
plt.fill_between(temp_time, master_contin5_flux[2238:]-master_contin5_error[2238:]-2.54e-12,master_contin5_flux[2238:]+ master_contin5_error[2238:]-2.54e-12, color = 'magenta',alpha = 0.25)

# plt.plot(master_time[2238:], master_He2_flux[2238:], color = 'red', linewidth = 1)
# plt.fill_between(master_time[2238:], master_He2_flux[2238:]-master_He2_error[2238:],master_He2_flux[2238:]+ master_He2_error[2238:], color = 'red',alpha = 0.25)
# plt.ylim(-2.5e-12,5e-12)
plt.legend()
plt.xlim(8500, 10500) #master_time[2658:2718]
# plt.ylim(.5e-12, 2.2e-12)
plt.xlabel(r'Time (s)', fontsize = 18)
plt.ylabel(r'Flux  ($\dfrac{erg}{cm^2s}$)', fontsize = 18)
plt.title("April 2021 Obs of TW Hya")





temp_time = np.linspace(1, len(master_time[2238:]),len(master_time[2238:]))*20
filename_counter = 72

plt.figure()

temp_change_in_obs_flag = (change_in_obs_flag - 44.76)*1000
for i in temp_change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > -1:
        plt.axvline(x=i, color = 'black', alpha = 0.5)
        # plt.text(x=i,y=-2.3e-12,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        # plt.text(x=i+.13,y=-2.3e-12,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        # filename_counter+=1

# plt.plot(temp_time, master_si4_flux[2238:], color = 'green', linewidth = 1)
# plt.fill_between(temp_time, master_si4_flux[2238:]-master_si4_error[2238:],master_si4_flux[2238:]+ master_si4_error[2238:], color = 'green',alpha = 0.25)
plt.plot(temp_time, master_c4_flux[2238:]/master_c4_flux[2690], color = 'green', linewidth = 1, label = 'C IV')
# plt.fill_between(temp_time, master_c4_flux[2238:]-master_c4_error[2238:]-3.28e-12,master_c4_flux[2238:]+ master_c4_error[2238:]-3.28e-12,color = 'green', alpha = 0.25)

plt.plot(temp_time, master_H2_flux[2238:]/master_H2_flux[2690], color = 'red', linewidth = 1, label = 'H2')
# plt.fill_between(temp_time, master_H2_flux[2238:]-master_H2_error[2238:]-1.23e-12,master_H2_flux[2238:]+ master_H2_error[2238:]-1.23e-12, color = 'red',alpha = 0.25)

plt.plot(temp_time, master_contin3_flux[2238:]/master_contin3_flux[2690], color = 'blue', linewidth = 1, label = 'Continuum in (1430 - 1520)$\AA$')
# plt.fill_between(temp_time, master_contin3_flux[2238:]-master_contin3_error[2238:]-2.06e-12,master_contin3_flux[2238:]+ master_contin3_error[2238:]-2.06e-12, color = 'blue',alpha = 0.25)

plt.plot(temp_time, master_contin4_flux[2238:]/master_contin4_flux[2690], color = 'orange', linewidth = 1, label = 'Continuum in (1600 - 1655)$\AA$')
# plt.fill_between(temp_time, master_contin4_flux[2238:]-master_contin4_error[2238:]+.43e-12,master_contin4_flux[2238:]+ master_contin4_error[2238:]+.43e-12, color = 'orange',alpha = 0.25)

plt.plot(temp_time, master_contin5_flux[2238:]/master_contin5_flux[2690], color = 'magenta', linewidth = 1, label = 'Continuum in (1680 - 1750)$\AA$')
# plt.fill_between(temp_time, master_contin5_flux[2238:]-master_contin5_error[2238:]-2.54e-12,master_contin5_flux[2238:]+ master_contin5_error[2238:]-2.54e-12, color = 'magenta',alpha = 0.25)

# plt.plot(master_time[2238:], master_He2_flux[2238:], color = 'red', linewidth = 1)
# plt.fill_between(master_time[2238:], master_He2_flux[2238:]-master_He2_error[2238:],master_He2_flux[2238:]+ master_He2_error[2238:], color = 'red',alpha = 0.25)
# plt.ylim(-2.5e-12,5e-12)
plt.legend()
plt.xlim(9060, 9350) #master_time[2658:2718]
plt.ylim(.8, 1.4)
plt.xlabel(r'Time (s)', fontsize = 18)
plt.ylabel(r'Flux  ($\dfrac{erg}{cm^2s}$)', fontsize = 18)
plt.title("Curves divided by flux(t=9060)")


#just looking at new observations
plt.figure()
plt.plot(temp_time, master_si4_flux[2238:], color = 'green', linewidth = 1)
plt.fill_between(temp_time, master_si4_flux[2238:]-master_si4_error[2238:],master_si4_flux[2238:]+ master_si4_error[2238:], color = 'green',alpha = 0.25, label = 'S IV')
plt.plot(temp_time, master_c4_flux[2238:], color = 'blue', linewidth = 1)
plt.fill_between(temp_time, master_c4_flux[2238:]-master_c4_error[2238:],master_c4_flux[2238:]+ master_c4_error[2238:],color = 'blue', alpha = 0.25, label = 'C IV')
plt.plot(temp_time, master_He2_flux[2238:], color = 'red', linewidth = 1)
plt.fill_between(temp_time, master_He2_flux[2238:]-master_He2_error[2238:],master_He2_flux[2238:]+ master_He2_error[2238:], color = 'red',alpha = 0.25, label = 'He II')
# plt.ylim(-2.5e-12,5e-12)
plt.xlabel(r'Time (ks)', fontsize = 18)
plt.ylabel(r'Flux  ($\dfrac{erg}{cm^2s}$)', fontsize = 18)
plt.title("April 2021 Obs of TW Hya")


























#anlayzinf commonBP




#mixing all of the ion species in one plot
fig,subplot = plt.subplots(3, 1, figsize = [19.3,9.5])
filename_counter = 0

for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i < master_time[1116] and flag == False:
        subplot[0].axvline(x=i, color = 'black', alpha = 0.5)
        # subplot[0].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        # subplot[0].text(x=i+.23,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i < master_time[1116] and flag == True:
        subplot[0].axvline(x=i, color = 'black', alpha = 0.5)
        # subplot[0].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        # subplot[0].text(x=i+.23,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1

# subplot[0].plot(master_time[0:1116], 1e12*master_si4_flux[0:1116], color = 'green', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
# subplot[0].fill_between(master_time[0:1116], 1e12*master_si4_flux[0:1116]-1e12*master_si4_error[0:1116],1e12*master_si4_flux[0:1116]+ 1e12*master_si4_error[0:1116], color = 'green',alpha = 0.25)    
# subplot[0].plot(master_time[0:1116], 1e12*master_c4_flux[0:1116], color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
# subplot[0].fill_between(master_time[0:1116], 1e12*master_c4_flux[0:1116]-1e12*master_c4_error[0:1116],1e12*master_c4_flux[0:1116]+ 1e12*master_c4_error[0:1116],color = 'blue',alpha = 0.25)
# subplot[0].plot(master_time[0:1116], 1e12*master_He2_flux[0:1116], color = 'red', linewidth = 1,label = r'He II 1640$\AA$')
# subplot[0].fill_between(master_time[0:1116], 1e12*master_He2_flux[0:1116]-1e12*master_He2_error[0:1116],1e12*master_He2_flux[0:1116]+ 1e12*master_He2_error[0:1116], color = 'red',alpha = 0.25)
subplot[0].plot(master_time[0:1116], 1e12*master_commonBP_flux[0:1116], color = 'black', linewidth = 1,label = r'CBP (1439 <$\lambda$ < 1561, 1635 <$\lambda$ < 1700) $\AA$')
subplot[0].fill_between(master_time[0:1116], 1e12*master_commonBP_flux[0:1116]-1e12*master_commonBP_error[0:1116],1e12*master_commonBP_flux[0:1116]+ 1e12*master_commonBP_error[0:1116], color = 'black',alpha = 0.25)

leg = subplot[0].legend(fontsize = 12, loc = 'upper right', markerscale =12)
# leg.get_lines()[0].set_linewidth(6)
# leg.get_lines()[1].set_linewidth(6)
# leg.get_lines()[2].set_linewidth(6)
# subplot[0].set_ylim(1e12*(-2.5e-12),1e12*5e-12)
subplot[0].set_xlim(-.25,22.5)

# labels = [item.get_text() for item in subplot[0].get_yticklabels()]
# labels[0] = ''; labels[1] = '';
# # labels[1] = 0;labels[2] = 1;labels[3] = 2;labels[4] = 3;labels[5] = 4;
# subplot[0].set_yticklabels(labels)


for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > master_time[1115] and i < master_time[2238] and flag == False:
        subplot[1].axvline(x=i, color = 'black', alpha = 0.5)
        # subplot[1].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        # subplot[1].text(x=i+.23,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i > master_time[1115] and i < master_time[2238] and flag == True:
        subplot[1].axvline(x=i, color = 'black', alpha = 0.5)
        # subplot[1].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        # subplot[1].text(x=i+.23,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1

# subplot[1].plot(master_time[1116:2238], 1e12*master_si4_flux[1116:2238], color = 'green', linewidth = 1)
# subplot[1].fill_between(master_time[1116:2238], 1e12*master_si4_flux[1116:2238]-1e12*master_si4_error[1116:2238],1e12*master_si4_flux[1116:2238]+ 1e12*master_si4_error[1116:2238], color = 'green',alpha = 0.25)
# subplot[1].plot(master_time[1116:2238], 1e12*master_c4_flux[1116:2238], color = 'blue', linewidth = 1)
# subplot[1].fill_between(master_time[1116:2238], 1e12*master_c4_flux[1116:2238]-1e12*master_c4_error[1116:2238],1e12*master_c4_flux[1116:2238]+ 1e12*master_c4_error[1116:2238],color = 'blue', alpha = 0.25)
# subplot[1].plot(master_time[1116:2238], 1e12*master_He2_flux[1116:2238], color = 'red', linewidth = 1)
# subplot[1].fill_between(master_time[1116:2238], 1e12*master_He2_flux[1116:2238]-1e12*master_He2_error[1116:2238],1e12*master_He2_flux[1116:2238]+ 1e12*master_He2_error[1116:2238], color = 'red',alpha = 0.25)
subplot[1].plot(master_time[1116:2238], 1e12*master_commonBP_flux[1116:2238], color = 'black', linewidth = 1)
subplot[1].fill_between(master_time[1116:2238], 1e12*master_commonBP_flux[1116:2238]-1e12*master_commonBP_error[1116:2238],1e12*master_commonBP_flux[1116:2238]+ 1e12*master_commonBP_error[1116:2238], color = 'black',alpha = 0.25)

# subplot[1].set_ylim(1e12*(-2.5e-12),1e12*5e-12)
subplot[1].set_xlim(22.1,45)

# labels = [item.get_text() for item in subplot[1].get_yticklabels()]
# labels[1] = 0;labels[2] = 1;labels[3] = 2;labels[4] = 3;labels[5] = 4;
# subplot[1].set_yticklabels(labels)


for i in change_in_obs_flag:
    flag = i in subplot0_HST_ID_fixes_list
    if i > master_time[2237] and flag == False:
        subplot[2].axvline(x=i, color = 'black', alpha = 0.5)
        # subplot[2].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        # subplot[2].text(x=i+.13,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1
    if i > master_time[2237] and flag == True:
        subplot[2].axvline(x=i, color = 'black', alpha = 0.5)
        # subplot[2].text(x=i,y=-2.3,s=folder_list[filename_counter][0:9], color = 'black', alpha = 0.75, rotation = 'vertical')
        # subplot[2].text(x=i+.13,y=-2.3,s=str(round(MJD_list[filename_counter],2)), color = 'black', alpha = 0.75, rotation = 'vertical')
        filename_counter+=1

# subplot[2].plot(master_time[2238:], 1e12*master_si4_flux[2238:], color = 'green', linewidth = 1)
# subplot[2].fill_between(master_time[2238:], 1e12*master_si4_flux[2238:]-1e12*master_si4_error[2238:],1e12*master_si4_flux[2238:]+ 1e12*master_si4_error[2238:], color = 'green',alpha = 0.25)
# subplot[2].plot(master_time[2238:], 1e12*master_c4_flux[2238:], color = 'blue', linewidth = 1)
# subplot[2].fill_between(master_time[2238:], 1e12*master_c4_flux[2238:]-1e12*master_c4_error[2238:],1e12*master_c4_flux[2238:]+ 1e12*master_c4_error[2238:],color = 'blue', alpha = 0.25)
# subplot[2].plot(master_time[2238:], 1e12*master_He2_flux[2238:], color = 'red', linewidth = 1)
# subplot[2].fill_between(master_time[2238:], 1e12*master_He2_flux[2238:]-1e12*master_He2_error[2238:],1e12*master_He2_flux[2238:]+ 1e12*master_He2_error[2238:], color = 'red',alpha = 0.25)
subplot[2].plot(master_time[2238:], 1e12*master_commonBP_flux[2238:], color = 'black', linewidth = 1)
subplot[2].fill_between(master_time[2238:], 1e12*master_commonBP_flux[2238:]-1e12*master_commonBP_error[2238:],1e12*master_commonBP_flux[2238:]+ 1e12*master_commonBP_error[2238:], color = 'black',alpha = 0.25)

# subplot[2].set_ylim(1e12*(-2.5e-12),1e12*5e-12)
subplot[2].set_xlim(44.65,56.85)


fig.text(0.5, 0.04, 'Time (ks)', ha='center', va='center', fontsize = 17)
fig.text(0.095, 0.5, r'Flux  ($erg \ cm^{-2} s^{-1}*10^{-12}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
fig.text(0.5, 0.92, 'COS G160M Observations of TW Hydrae ($\Delta t = %ss$, 15.75hrs Total)' % time_resolution, ha='center', va='center', fontsize = 18)





