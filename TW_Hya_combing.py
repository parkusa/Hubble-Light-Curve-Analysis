
# -*- coding: utf-8 -*-
"""
Taking the TW data and combing it for possible flare events
"""
import math
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from astropy.io import fits
from tqdm import tqdm
import time as ostime
from scipy.optimize import curve_fit


def Reverse(lst):
    return [ele for ele in reversed(lst)]

def find_nearest(array,value):
    idx = np.nanargmin(np.abs(array - value))
    return idx


def subtract_quies(species_count_rate, quiescent_time_index):
    mean_quies = np.mean(species_count_rate[quiescent_time_index])
    subtracted_counts = species_count_rate - mean_quies
    return subtracted_counts

def onclick(event):
    global ix, iy, coords_x, coords_y, terminator
    ix, iy = event.xdata, event.ydata
    coords_x.append(ix)
    coords_y.append(iy)
    print('time:  '+ str(ix))
    if len(coords_x) == 1:
        fig.canvas.mpl_disconnect(fig)
        plt.close()
        coords_x = np.sort(coords_x)
        terminator = 1
    return
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

def func_linear(x, m, b):
    return (m*x+b)

#defining the time resolution of this file
time_resolution = 20
#getting the name of all folders in /out, each one represents a full observation,
#the files inside are subframes
folder_list = get_file_names_with_strings(['l']) #all of the files
# This block of code was used to verify that the files I read in are in MJD
# sorted order
MJD_list = [] #defining empty MJD list, to be filled in 'for loop'
for i in folder_list:
    filename = i + '-total-x1d.fits'
    # print(filename)
    hdul = fits.open(r"C:\Windows\System32\astroconda\astroconda_TW_Hya\COS_data\out20\%s\%s" % (i,filename))
    TEXPSTRT = hdul[1].header['EXPSTART'] #start time (MJD) of 1st exposure in file
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
master_commonBP_flux=np.asarray([]); master_commonBP_counts=np.asarray([]); master_commonBP_error=np.asarray([]);


master_time=[]; time = 0
change_in_obs_flag = [] 


folder_count = 1
for i in tqdm(folder_list): # this loop iterates through each folder, a nested loop
#iterates through subframes and generates lightcurve
    print(' ');print(' ');
    print('Now working on %s of 72 G160 folders' % folder_count)
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
            ##############################################################
            
            #Now to compute the average fluxes and count rates from total x1d file#
    
            #avg flux
            avg_c2_line_flux = np.trapz(flux[c2_line_location],wave[c2_line_location])
            avg_c3_line_flux = np.trapz(flux[c3_line_location],wave[c3_line_location])
            avg_n5_line_flux = np.trapz(flux[n5_line_location],wave[n5_line_location])
            avg_si4_line_flux = np.trapz(flux[si4_line_location],wave[si4_line_location])
            avg_c4_line_flux = np.trapz(flux[c4_line_location],wave[c4_line_location])
            avg_He2_line_flux = np.trapz(flux[He2_line_location],wave[He2_line_location])
            avg_commonBP_line_flux = np.trapz(flux[commonBP_line_location],wave[commonBP_line_location])
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
            avg_commonBP_line_counts = np.sum(counts[commonBP_line_location])
            avg_continuum_line_counts = np.sum(counts[continuum_line_location])
            avg_bandpass1_line_counts = np.sum(counts[bandpass1_line_location])
            avg_bandpass2_line_counts = np.sum(counts[bandpass2_line_location])
            
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
        commonBP_flux=[]; commonBP_flux_error=[]; commonBP_counts=[]; commonBP_error=[]
        #adding extra lines to compare with Hawley et al 2003
        c2_flux=[]; c2_flux_error=[]; c2_counts=[]; c2_error=[];
        si2_flux=[]; si2_flux_error=[]; si2_counts=[]; si2_error=[];
        si3_flux=[]; si3_flux_error=[]; si3_counts=[]; si3_error=[];
        o1_flux=[]; o1_flux_error=[]; o1_counts=[]; o1_error=[];
        #larger non-line regions
        continuum_flux=[]; continuum_flux_error=[]; continuum_counts=[]; continuum_error=[];
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
        ################################################################
        
        # computing count rates
        c2temp_counts = np.sum(counts[c2_line_location])
        c3temp_counts = np.sum(counts[c3_line_location])
        c4temp_counts = np.sum(counts[c4_line_location])
        
        continuum_temp_counts = np.sum(counts[continuum_line_location])
        si4_temp_counts = np.sum(counts[si4_line_location])
        He2_temp_counts = np.sum(counts[He2_line_location])
        commonBP_temp_counts = np.sum(counts[commonBP_line_location])
        
        #appending count rates
        c2_counts.append(c2temp_counts); c2_error.append(np.sqrt(c2temp_counts*time_resolution)/time_resolution)
        c3_counts.append(c3temp_counts); c3_error.append(np.sqrt(c3temp_counts*time_resolution)/time_resolution)
        c4_counts.append(c4temp_counts); c4_error.append(np.sqrt(c4temp_counts*time_resolution)/time_resolution)
        
        continuum_counts.append(continuum_temp_counts); continuum_error.append(np.sqrt(continuum_temp_counts*time_resolution)/time_resolution)
        si4_counts.append(si4_temp_counts); si4_error.append(np.sqrt(si4_temp_counts*time_resolution)/time_resolution)
        He2_counts.append(He2_temp_counts); He2_error.append(np.sqrt(He2_temp_counts*time_resolution)/time_resolution)
        commonBP_counts.append(commonBP_temp_counts); commonBP_error.append(np.sqrt(commonBP_temp_counts*time_resolution)/time_resolution)
        
        
        #converting count lists to arrays
        c2_counts=np.asarray(c2_counts); c2_error=np.asarray(c2_error);
        c3_counts=np.asarray(c3_counts); c3_error=np.asarray(c3_error);
        c4_counts=np.asarray(c4_counts); c4_error=np.asarray(c4_error);
        
        continuum_counts=np.asarray(continuum_counts); continuum_error=np.asarray(continuum_error);
        si4_counts=np.asarray(si4_counts); si4_error=np.asarray(si4_error);
        He2_counts=np.asarray(He2_counts); He2_error=np.asarray(He2_error);
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
del c3_line_location;del c3temp_counts;

#scaling
master_time = np.asarray(master_time); change_in_obs_flag = np.asarray(change_in_obs_flag);
master_time = master_time/1000; change_in_obs_flag = change_in_obs_flag/1000

change_in_obs_index = []
for i in change_in_obs_flag:
    temp = np.where(master_time == i)
    temp = temp[0][0]
    change_in_obs_index.append(temp)
    del temp





# target_func = func_linear
# three_sigma_list = []
# three_sigma_list_err_up = []
# three_sigma_list_err_down = []
# five_sigma_list_err_up = []
# five_sigma_list_err_down = []
# mean_flux_list = []
# std_deviation_list = []
# time_observed = 0
# for i in range(len(change_in_obs_index)):
#     circuit_breaker =1
#     while circuit_breaker ==1:
#         if i == np.max(range(len(change_in_obs_index))):
#             time = master_time[change_in_obs_index[i]:-1]
#         else:
#             time = master_time[change_in_obs_index[i]:change_in_obs_index[i+1]]
        
        
#         ###   pick which line you're analyzing ###
#         #________________________________________#
        
#         if i == np.max(range(len(change_in_obs_index))):
#             flux = master_si4_flux[change_in_obs_index[i]:-1]
#             err  = master_si4_error[change_in_obs_index[i]:-1]
#         else:
#             flux = master_si4_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
#             err  = master_si4_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
#         species_name = 'S IV'
        
#         # if i == np.max(range(len(change_in_obs_index))):
#         #     flux = master_c4_flux[change_in_obs_index[i]:-1]
#         #     err  = master_c4_error[change_in_obs_index[i]:-1]
#         # else:
#         #     flux = master_c4_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
#         #     err  = master_c4_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
#         # species_name = 'C IV'
        
#         # if i == np.max(range(len(change_in_obs_index))):
#         #     flux = master_He2_flux[change_in_obs_index[i]:-1]
#         #     err  = master_He2_error[change_in_obs_index[i]:-1]
#         # else:
#         #     flux = master_He2_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
#         #     err  = master_He2_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
#         # species_name = 'He II'
        

#         # if len(flux) < 4: break
#         # if type(flux[0]) !=
#        # ############ This section finds all data above 3 and 5 sigma levels
#         if math.isnan(flux[0]): break
    
#         mean_flux_list.append(np.mean(flux))
#         std_deviation_list.append(np.std(flux))
#         time_observed+=20*len(time)


       
#         slope_guess = (flux[-1] - flux[0])/(len(time)*(time[1]-time[0]))
#         popt, pcov = curve_fit(target_func, time, flux, maxfev=2000,
#                                 p0 = np.asarray([slope_guess,0]))
#         fit_line = target_func(np.array(time), *popt)
        
#         adjusted_flux = (np.array(range(len(flux)))*(-1)*popt[0]*(time[1]-time[0])
#                       + flux)
        
        
#         # don't consider data points > 1.7*median when finding 3 and 5 sigma levels
#         # and their two neighboring points
#         temp_flux = adjusted_flux
#         temp_time = time
#         #'''
#         temp_potential_flare = list(np.where(temp_flux > 1.7*np.median(temp_flux))[0])
#         potential_flare = []
#         for j in temp_potential_flare:
#             if j > 0:
#                 potential_flare.insert(0, j-1)
#             potential_flare.insert(0, j)
#             if j < len(flux):
#                 potential_flare.insert(0, j+1)
#         potential_flare = list(dict.fromkeys(potential_flare)) # converting to a dictionary gets rid of
#                                                               #duplicates. Then back to a list for uses
#         for j in potential_flare:
#             temp_flux = np.delete(temp_flux,j) 
#             temp_time = np.delete(temp_time,j)
            
            
#         # Now getting rid of all points > 3 sigma
#         temp_std = np.std(temp_flux)
#         temp_mean = np.mean(temp_flux)
#         temp_potential_flare = list(np.where(temp_flux > temp_mean + 2*temp_std)[0])
#         # temp_potential_flare = list(np.where(temp_flux > temp_mean + 2*temp_std)[0])
#         potential_flare = []
#         for j in temp_potential_flare:
#             if j > 0:
#                 potential_flare.insert(0, j-1)
#             potential_flare.insert(0, j)
#             if j < len(flux) - 1:
#                 potential_flare.insert(0, j+1)
#         potential_flare = list(dict.fromkeys(potential_flare)) # converting to a dictionary gets rid of
#                                                               #duplicates. Then back to a list for uses
#         for j in potential_flare:
#             temp_flux = np.delete(temp_flux,j) 
#             temp_time = np.delete(temp_time,j)
        
#     #'''
    
#         mean1 = np.mean(flux)
#         std1 = np.std(flux)
#         # mean = np.mean(adjusted_flux)
#         # std = np.std(adjusted_flux)
#         mean = np.mean(temp_flux)
#         std = np.std(temp_flux)
#         fig,subplot = plt.subplots(2, 1)
        
#         subplot[0].plot(time*1000, flux, marker = '.')
#         subplot[0].plot(time*1000, fit_line, color = 'red', marker = '.')
#         subplot[0].plot(time*1000, adjusted_flux, color = 'magenta')
#         subplot[0].plot(time*1000, np.ones(len(time*1000))*mean1+3*std1, color = 'orange', linestyle = 'dashed')
#         subplot[0].plot(time*1000, np.ones(len(time*1000))*mean1+5*std1, color = 'red')
        
#         subplot[1].errorbar(time*1000, adjusted_flux, yerr=err, color = 'magenta')
#         subplot[1].plot(time*1000, np.ones(len(time*1000))*mean, color = 'black')
#         subplot[1].plot(time*1000, np.ones(len(time*1000))*mean+3*std, color = 'orange', linestyle = 'dashed')
#         subplot[1].plot(time*1000, np.ones(len(time*1000))*mean+5*std, color = 'red')
        
#         fig.text(0.5, 0.04, 'Time (s)', ha='center', va='center', fontsize = 17)
#         fig.text(0.035, 0.5, r'$Flux$  ($erg*cm^{-2}s^{-1}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
#         fig.text(0.5, 0.92, 'TW Hya, %s, %s' % (folder_list[i],species_name), ha='center', va='center', fontsize = 18)
        
#         if i<10: im_name ='00%s' % i
#         if i>9 and i<100: im_name ='0%s' % i
#         if i>99: im_name ='%s' % i
        
#         # fig.savefig("im%s.jpg" % im_name) #save for movies :)
#         temp = np.where(adjusted_flux > mean+3*std)
#         temp = temp[0]
        
#         for j in temp: three_sigma_list.append(time[int(j)])
#         temp = np.where(adjusted_flux+err > mean+3*std)
#         temp = temp[0]
#         for j in temp: three_sigma_list_err_up.append(time[int(j)])
#         temp = np.where(adjusted_flux-err > mean+3*std)
#         temp = temp[0]
#         if len(temp) > 0: print('3 sigma flare detected in %s ! img # %s' % (folder_list[i], i))
#         for j in temp: three_sigma_list_err_down.append(time[int(j)])
#         temp = np.where(adjusted_flux+err > mean+5*std)
#         temp = temp[0]
#         if len(temp) > 0: print('potential 5 sigma flare detected in %s img # %s' % (folder_list[i], i))
#         for j in temp: five_sigma_list_err_up.append(time[int(j)])
#         temp = np.where(adjusted_flux-err > mean+5*std)
#         temp = temp[0]
#         if len(temp) > 0:
#             print('5 sigma flare detected in %s !!! img # %s' % (folder_list[i], i))
#             for j in temp:
#                 print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! %s" % im_name)
#                 five_sigma_list_err_down.append(time[int(j)])
#         ####### Uncomment the following code to see each plot individually
#         ####### and have to click to close it
#         '''
#         #doesn't actually work with subplots yet...'
        
#         coords_x = []; coords_y = [];
#         fig.canvas.mpl_connect('button_press_event', onclick)
#         print('click anywhere to close plot')
#         terminator = 0
#         while terminator == 0:
#             plt.pause(.001)
#         del coords_x;del coords_y;
#         '''
#         #######
#         #######
#         plt.close('all')
#         del temp;  del mean; del std;del adjusted_flux;
#         del fit_line; del slope_guess; del popt; del pcov;
#         del temp_flux; del temp_mean; del temp_potential_flare;
#         del temp_std; del temp_time;

#         circuit_breaker = 0


# three_sigma_index = []
# for i in three_sigma_list: three_sigma_index.append(np.where(master_time == i))
# three_sigma_index_err_up = []
# for i in three_sigma_list_err_up:three_sigma_index_err_up.append(np.where(master_time == i)[0][0])
# three_sigma_index_err_down = []
# for i in three_sigma_list_err_down:three_sigma_index_err_down.append(np.where(master_time == i)[0][0])
# five_sigma_index_err_up = []
# for i in five_sigma_list_err_up:five_sigma_index_err_up.append(np.where(master_time == i)[0][0])
# five_sigma_index_err_down = []
# for i in five_sigma_list_err_down:five_sigma_index_err_down.append(np.where(master_time == i)[0][0])



# print('\n')
# print('For %s:' % species_name)
# temp1 = (len(five_sigma_list_err_up)-len(five_sigma_list_err_down))/2 +len(five_sigma_list_err_down)
# temp2 = (len(five_sigma_list_err_up)-len(five_sigma_list_err_down))/2
# print(">5 sig : %s +/- %s, which goes from %s to %s" % (temp1,temp2,temp1 - temp2, temp1 + temp2))
# temp1 = (len(three_sigma_list_err_up)-len(three_sigma_list_err_down))/2 +len(three_sigma_list_err_down)
# temp2 = (len(three_sigma_list_err_up)-len(three_sigma_list_err_down))/2
# print(">3 sig : %s +/- %s, which goes from %s to %s" % (temp1,temp2,temp1 - temp2, temp1 + temp2))

# print('Time observed = %s s' % time_observed)






# np.savetxt('TW_Hya_flux_stats.csv', np.column_stack((mean_flux_list, std_deviation_list)), delimiter = ',', fmt = '%s')




print('\n')
print("now computing the estimates using the Hilton 2011 approach:")
print('\n')


def Hilton_LCA(adjusted_flux, mean, std, time,i,announce):
    three_sigma_list_temp = []
    five_sigma_list_temp = []
    observed_data = adjusted_flux
    sigma_list = (observed_data-mean)/std
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
                    print('\n')
                    # print('a3 sigma flare detected in %s at %s ! img # %s' % (folder_list[i],time_of_flare_peak, i))
        #deleting flares which don't breach the 5sigma level
        flares_list[:] = [sublist for sublist in flares_list if any(5 <= sigma_list[x] for x in sublist)]            
        #appending and decaring 5 sigma flares
        if len(flares_list) > 0:                
            for x in flares_list:
                loc_of_flare_peak = np.where(sigma_list == np.max(sigma_list[x]))[0]
                time_of_flare_peak = time[loc_of_flare_peak]
                five_sigma_list_temp.append(time_of_flare_peak)
                if announce == 'on':
                    print('b5 sigma flare detected in %s at %s !!!!! img # %s' % (folder_list[i],time_of_flare_peak, i))
                
    if num_of_individ_flares == 'list_of_one' and len(flares_list[0]) >= 3:
        max_sig = np.max(sigma_list[tuple(flares_list)])
        if max_sig > 3:                               
            loc_of_flare_peak = np.where(sigma_list == max_sig)[0]
            time_of_flare_peak = time[loc_of_flare_peak]
            three_sigma_list_temp.append(time_of_flare_peak)
            if announce == 'on':
                print('\n')
                # print('c3 sigma flare detected in %s at %s ! img # %s' % (folder_list[i],time_of_flare_peak, i))
        if max_sig > 5:
            loc_of_flare_peak = np.where(sigma_list == max_sig)[0]
            time_of_flare_peak = time[loc_of_flare_peak]
            five_sigma_list_temp.append(time_of_flare_peak)
            if announce == 'on':
                print('d5 sigma flare detected in %s at %s !!!!! img # %s' % (folder_list[i],time_of_flare_peak, i))



    return three_sigma_list_temp, five_sigma_list_temp






target_func = func_linear
three_sigma_list = []
three_sigma_list_err_up = []
three_sigma_list_err_down = []
five_sigma_list = []
five_sigma_list_err_up = []
five_sigma_list_err_down = []
mean_flux_list = []
std_deviation_list = []
time_observed = 0
for i in range(len(change_in_obs_index)):
    circuit_breaker =1
    while circuit_breaker ==1:
        if i == np.max(range(len(change_in_obs_index))):
            time = master_time[change_in_obs_index[i]:-1]
        else:
            time = master_time[change_in_obs_index[i]:change_in_obs_index[i+1]]
        
        
        ###   pick which line you're analyzing ###
        #________________________________________#
        
        if i == np.max(range(len(change_in_obs_index))):
            flux = master_commonBP_flux[change_in_obs_index[i]:-1]
            err  = master_commonBP_error[change_in_obs_index[i]:-1]
        else:
            flux = master_commonBP_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
            err  = master_commonBP_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
        species_name = 'commonBP'
        
        
        # if i == np.max(range(len(change_in_obs_index))):
        #     flux = master_si4_flux[change_in_obs_index[i]:-1]
        #     err  = master_si4_error[change_in_obs_index[i]:-1]
        # else:
        #     flux = master_si4_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
        #     err  = master_si4_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
        # species_name = 'S IV'
        
        # if i == np.max(range(len(change_in_obs_index))):
        #     flux = master_c4_flux[change_in_obs_index[i]:-1]
        #     err  = master_c4_error[change_in_obs_index[i]:-1]
        # else:
        #     flux = master_c4_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
        #     err  = master_c4_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
        # species_name = 'C IV'
        
        # if i == np.max(range(len(change_in_obs_index))):
        #     flux = master_He2_flux[change_in_obs_index[i]:-1]
        #     err  = master_He2_error[change_in_obs_index[i]:-1]
        # else:
        #     flux = master_He2_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
        #     err  = master_He2_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
        # species_name = 'He II'
        

        # if len(flux) < 4: break
        # if type(flux[0]) !=
       # ############ This section finds all data above 3 and 5 sigma levels
        if math.isnan(flux[0]): break
    
        mean_flux_list.append(np.mean(flux))
        std_deviation_list.append(np.std(flux))
        time_observed+=20*len(time)


       
        slope_guess = (flux[-1] - flux[0])/(len(time)*(time[1]-time[0]))
        popt, pcov = curve_fit(target_func, time, flux, maxfev=2000,
                                p0 = np.asarray([slope_guess,0]))
        fit_line = target_func(np.array(time), *popt)
        
        adjusted_flux = (np.array(range(len(flux)))*(-1)*popt[0]*(time[1]-time[0])
                      + flux)
        
        
        # don't consider data points > 1.7*median when finding 3 and 5 sigma levels
        # and their two neighboring points
        temp_flux = adjusted_flux
        temp_time = time
        #'''
        temp_potential_flare = list(np.where(temp_flux > 1.7*np.median(temp_flux))[0])
        potential_flare = []
        for j in temp_potential_flare:
            if j > 0:
                potential_flare.insert(0, j-1)
            potential_flare.insert(0, j)
            if j < len(flux):
                potential_flare.insert(0, j+1)
        potential_flare = list(dict.fromkeys(potential_flare)) # converting to a dictionary gets rid of
                                                              #duplicates. Then back to a list for uses
        for j in potential_flare:
            temp_flux = np.delete(temp_flux,j) 
            temp_time = np.delete(temp_time,j)
            
            
        # Now getting rid of all points > 3 sigma
        temp_std = np.std(temp_flux)
        temp_mean = np.mean(temp_flux)
        temp_potential_flare = list(np.where(temp_flux > temp_mean + 2*temp_std)[0])
        # temp_potential_flare = list(np.where(temp_flux > temp_mean + 2*temp_std)[0])
        potential_flare = []
        for j in temp_potential_flare:
            if j > 0:
                potential_flare.insert(0, j-1)
            potential_flare.insert(0, j)
            if j < len(flux) - 1:
                potential_flare.insert(0, j+1)
        potential_flare = list(dict.fromkeys(potential_flare)) # converting to a dictionary gets rid of
                                                              #duplicates. Then back to a list for uses
        for j in potential_flare:
            temp_flux = np.delete(temp_flux,j) 
            temp_time = np.delete(temp_time,j)
        
    #'''
    
        mean1 = np.mean(flux)
        std1 = np.std(flux)
        # mean = np.mean(adjusted_flux)
        # std = np.std(adjusted_flux)
        mean = np.mean(temp_flux)
        std = np.mean(err)
        fig,subplot = plt.subplots(2, 1)
        
        subplot[0].plot(time*1000, flux, marker = '.')
        subplot[0].plot(time*1000, fit_line, color = 'red', marker = '.')
        subplot[0].plot(time*1000, adjusted_flux, color = 'magenta')
        subplot[0].plot(time*1000, np.ones(len(time*1000))*mean1+3*std1, color = 'orange', linestyle = 'dashed')
        subplot[0].plot(time*1000, np.ones(len(time*1000))*mean1+5*std1, color = 'red')
        
        subplot[1].errorbar(time*1000, adjusted_flux, yerr=err, color = 'magenta')
        subplot[1].plot(time*1000, np.ones(len(time*1000))*mean, color = 'black')
        subplot[1].plot(time*1000, np.ones(len(time*1000))*mean+3*std, color = 'orange', linestyle = 'dashed')
        subplot[1].plot(time*1000, np.ones(len(time*1000))*mean+5*std, color = 'red')
        
        fig.text(0.5, 0.04, 'Time (s)', ha='center', va='center', fontsize = 17)
        fig.text(0.035, 0.5, r'$Flux$  ($erg*cm^{-2}s^{-1}$)', ha='center', va='center', rotation='vertical', fontsize = 17)
        fig.text(0.5, 0.92, 'TW Hya, %s, %s' % (folder_list[i],species_name), ha='center', va='center', fontsize = 18)
        
        if i<10: im_name ='00%s' % i
        if i>9 and i<100: im_name ='0%s' % i
        if i>99: im_name ='%s' % i
        
        
        std_of_each_datum = (adjusted_flux-mean)/std
        # fig.savefig("im%s.jpg" % im_name) #save for movies :)
        
        plt.close('all')

    
        # announce = 'on'
        # observed_data = adjusted_flux + err #change +,-, or no error depending on what you're checking
        # sigma_list = (observed_data-mean)/std
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
        #             time_of_flare_peak = time[loc_of_flare_peak]
        #             three_sigma_list.append(time_of_flare_peak)
        #             if announce == 'on':
        #                 print('a3 sigma flare detected in %s at %s ! img # %s' % (folder_list[i],time_of_flare_peak, i))
        #     #deleting flares which don't breach the 5sigma level
        #     flares_list[:] = [sublist for sublist in flares_list if any(5 <= sigma_list[x] for x in sublist)]            
        #     #appending and decaring 5 sigma flares
        #     if len(flares_list) > 0:                
        #         for x in flares_list:
        #             loc_of_flare_peak = np.where(sigma_list == np.max(sigma_list[x]))[0]
        #             time_of_flare_peak = time[loc_of_flare_peak]
        #             five_sigma_list.append(time_of_flare_peak)
        #             if announce == 'on':
        #                 print('b5 sigma flare detected in %s at %s !!!!! img # %s' % (folder_list[i],time_of_flare_peak, i))
                    
        # if num_of_individ_flares == 'list_of_one' and len(flares_list[0]) >= 3:
        #     max_sig = np.max(sigma_list[tuple(flares_list)])
        #     if max_sig > 3:                               
        #         loc_of_flare_peak = np.where(sigma_list == max_sig)[0]
        #         time_of_flare_peak = time[loc_of_flare_peak]
        #         three_sigma_list.append(time_of_flare_peak)
        #         if announce == 'on':
        #             print('c3 sigma flare detected in %s at %s ! img # %s' % (folder_list[i],time_of_flare_peak, i))
        #     if max_sig > 5:
        #         loc_of_flare_peak = np.where(sigma_list == max_sig)[0]
        #         time_of_flare_peak = time[loc_of_flare_peak]
        #         five_sigma_list.append(time_of_flare_peak)
        #         if announce == 'on':
        #             print('d5 sigma flare detected in %s at %s !!!!! img # %s' % (folder_list[i],time_of_flare_peak, i))


        a,b = Hilton_LCA(adjusted_flux, mean, std, time, i,announce = 'on')
        for x in a:
            three_sigma_list.append(x)
        for x in b:
            five_sigma_list.append(x)
       
        a,b = Hilton_LCA(adjusted_flux+err, mean, std, time, i,announce = 'off')
        for x in a:
            three_sigma_list_err_up.append(x)
        for x in b:
            five_sigma_list_err_up.append(x)            
        
        a,b = Hilton_LCA(adjusted_flux-err, mean, std, time, i,announce = 'on')
        for x in a:
            three_sigma_list_err_down.append(x)
        for x in b:
            five_sigma_list_err_down.append(x)            



        circuit_breaker = 0
print('\n')
print('For %s:' % species_name)
temp0 = len(five_sigma_list)
temp1 = len(five_sigma_list_err_up)-temp0
temp2 = temp0 - len(five_sigma_list_err_down)
print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
temp0 = len(three_sigma_list)
temp1 = len(three_sigma_list_err_up)-temp0
temp2 = temp0 - len(three_sigma_list_err_down)
print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))


print('Time observed = %s s' % time_observed)



        ######


               

        
        
        
        
        
        ### first attemot at re-creatign hilton 2011
        
#         observed_data = adjusted_flux
#         for x in range(len(observed_data)-2):
#             #finding base 3 sigma events
#             temp = observed_data[x]
#             if temp > 3*std+mean:
#                 temp2 = observed_data[x+1];temp3 = observed_data[x+2];
#                 if temp2 > 2.5*std+mean and temp3 > 2.5*std+mean:
#                     if time[x] not in three_sigma_list and time[x-1] not in three_sigma_list and time[x-2] not in three_sigma_list:
#                         three_sigma_list.append(time[x])
#             #finding base 5 sigma events
#             if temp > 5*std+mean:
#                 temp2 = observed_data[x+1];temp3 = observed_data[x+2];
#                 if temp2 > 2.5*std+mean and temp3 > 2.5*std+mean:
#                     if time[x] not in five_sigma_list and time[x-1] not in five_sigma_list and time[x-2] not in five_sigma_list:
#                         five_sigma_list.append(time[x])    
#             #finding upper 3 sigma events
#             temp_err = err[x]
#             temp = observed_data[x] + temp_err
#             if temp > 3*std+mean:
#                 temp_err2 = err[x+1]
#                 temp_err3 = err[x+2]        
#                 temp2 = observed_data[x+1] + temp_err2;
#                 temp3 = observed_data[x+2] + temp_err3;
#                 if temp2 > 2.5*std+mean and temp3 > 2.5*std+mean:
#                     if time[x] not in three_sigma_list_err_up and time[x-1] not in three_sigma_list_err_up and time[x-2] not in three_sigma_list_err_up:
#                         three_sigma_list_err_up.append(time[x])
#             #finding lower 3 sigma events            
#             temp_err = err[x]
#             temp = observed_data[x] - temp_err
#             if temp > 3*std+mean:
#                 temp_err2 = err[x+1]
#                 temp_err3 = err[x+2]        
#                 temp2 = observed_data[x+1] - temp_err2;
#                 temp3 = observed_data[x+2] - temp_err3;
#                 if temp2 > 2.5*std+mean and temp3 > 2.5*std+mean:
#                     if time[x] not in three_sigma_list_err_down and time[x-1] not in three_sigma_list_err_down and time[x-2] not in three_sigma_list_err_down:
#                         three_sigma_list_err_down.append(time[x])
#             #finding upper 5 sigma events
#             temp_err = err[x]
#             temp = observed_data[x] + temp_err
#             if temp > 5*std+mean:
#                 temp_err2 = err[x+1]
#                 temp_err3 = err[x+2]        
#                 temp2 = observed_data[x+1] + temp_err2;
#                 temp3 = observed_data[x+2] + temp_err3;
#                 if temp2 > 2.5*std+mean and temp3 > 2.5*std+mean:
#                     if time[x] not in five_sigma_list_err_up and time[x-1] not in five_sigma_list_err_up and time[x-2] not in five_sigma_list_err_up:
#                         five_sigma_list_err_up.append(time[x])
#             #finding lower 5 sigma events        
#             temp_err = err[x]
#             temp = observed_data[x] - temp_err
#             if temp > 5*std+mean:
#                 temp_err2 = err[x+1]
#                 temp_err3 = err[x+2]        
#                 temp2 = observed_data[x+1] - temp_err2;
#                 temp3 = observed_data[x+2] - temp_err3;
#                 if temp2 > 2.5*std+mean and temp3 > 2.5*std+mean:
#                     if time[x] not in five_sigma_list_err_down and time[x-1] not in five_sigma_list_err_down and time[x-2] not in five_sigma_list_err_down:
#                         five_sigma_list_err_down.append(time[x])
# print('\n')
# print('For %s:' % species_name)
# temp0 = len(five_sigma_list)
# temp1 = len(five_sigma_list_err_up)-temp0
# temp2 = temp0 - len(five_sigma_list_err_down)
# print(">5 sig : %s + %s - %s" % (temp0, temp1,temp2))
# temp0 = len(three_sigma_list)
# temp1 = len(three_sigma_list_err_up)-temp0
# temp2 = temp0 - len(three_sigma_list_err_down)
# print(">3 sig : %s + %s - %s" % (temp0, temp1,temp2))


# print('Time observed = %s s' % time_observed)















#removing the most deviant point till std asymptotes strategy

# target_func = func_linear
# three_sigma_list = []
# three_sigma_list_err_up = []
# three_sigma_list_err_down = []
# five_sigma_list_err_up = []
# five_sigma_list_err_down = []
# mean_flux_list = []
# std_deviation_list = []


# for i in range(30):

#     circuit_breaker =1
#     while circuit_breaker ==1:
#         if i == np.max(range(len(change_in_obs_index))):
#             time = master_time[change_in_obs_index[i]:-1]
#         else:
#             time = master_time[change_in_obs_index[i]:change_in_obs_index[i+1]]
        
        
#         ###   pick which line you're analyzing ###
#         #________________________________________#
        
#         # if i == np.max(range(len(change_in_obs_index))):
#         #     flux = master_si4_flux[change_in_obs_index[i]:-1]
#         #     err  = master_si4_error[change_in_obs_index[i]:-1]
#         # else:
#         #     flux = master_si4_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
#         #     err  = master_si4_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
        
        
#         if i == np.max(range(len(change_in_obs_index))):
#             flux = master_c4_flux[change_in_obs_index[i]:-1]
#             err  = master_c4_error[change_in_obs_index[i]:-1]
#         else:
#             flux = master_c4_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
#             err  = master_c4_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
        
        
#         # if i == np.max(range(len(change_in_obs_index))):
#         #     flux = master_He2_flux[change_in_obs_index[i]:-1]
#         #     err  = master_He2_error[change_in_obs_index[i]:-1]
#         # else:
#         #     flux = master_He2_flux[change_in_obs_index[i]:change_in_obs_index[i+1]]
#         #     err  = master_He2_error[change_in_obs_index[i]:change_in_obs_index[i+1]]
        
    
#         # if len(flux) < 4: break
#         # if type(flux[0]) !=
#        # ############ This section finds all data above 3 and 5 sigma levels
#         if math.isnan(flux[0]): break
    
#         mean_flux_list.append(np.mean(flux))
#         std_deviation_list.append(np.std(flux))
    
    
    
       
#         slope_guess = (flux[-1] - flux[0])/(len(time)*(time[1]-time[0]))
#         popt, pcov = curve_fit(target_func, time, flux, maxfev=2000,
#                                 p0 = np.asarray([slope_guess,0]))
#         fit_line = target_func(np.array(time), *popt)
        
#         adjusted_flux = (np.array(range(len(flux)))*(-1)*popt[0]*(time[1]-time[0])
#                       + flux)
        
        
    
#         # don't consider data points > 1.7*median when finding 3 and 5 sigma levels
#         # and their two neighboring points
#         temp_flux = adjusted_flux
#         temp_time = time
        
        
#         std_list = []
#         it_count = []
#         for i in range(len(temp_flux)):
#             std = np.std(temp_flux)
#             std_list.append(std)
#             it_count.append(i)
#             x = (temp_flux-np.mean(temp_flux))/std
#             most_deviant_point = list(np.where(x == np.max(x))[0])
#             for j in most_deviant_point:
#                 temp_flux = np.delete(temp_flux,j) 
#                 temp_time = np.delete(temp_time,j)
        

#         circuit_breaker = 0
#     del adjusted_flux;
#     del fit_line; del slope_guess; del popt; del pcov;
#     del temp_flux; 
#     del temp_time;
#     plt.figure()
#     plt.plot(it_count,std_list,'.')
                    
                
    
    
    
    
    
    
    
    
    





#plots n styff?


# plt.figure()
# plt.plot(master_time, master_si4_flux, color = 'blue', linewidth = 1,label = r'Si IV (1394 + 1403)$\AA$')
# plt.fill_between(master_time, master_si4_flux-master_si4_error,master_si4_flux+ master_si4_error,color = 'blue',alpha = 0.25)
# plt.plot(master_time[three_sigma_index_err_down],master_si4_flux[three_sigma_index_err_down],
#          marker = 'x',markersize=7, color = 'red', linewidth = 0)
# plt.plot(master_time[three_sigma_index_err_up],master_si4_flux[three_sigma_index_err_up],
#          marker = 'x',markersize=7, color = 'green', linewidth = 0)
# plt.xlabel('Time (ks)', fontsize = 17)
# plt.ylabel( r'Flux  ($\dfrac{erg}{cm^2s}$) *$10^{-12}$',fontsize = 17)
# plt.title('COS G160 Observations of Si IV TW Hydrae ($\Delta t = %ss$, 12.4hrs Total)' % time_resolution,fontsize = 18)




# plt.figure()
# plt.plot(master_time, master_c4_flux, color = 'blue', linewidth = 1,label = r'C IV (1548.2 + 1550.7)$\AA$')
# plt.fill_between(master_time, master_c4_flux-master_c4_error,master_c4_flux+ master_c4_error,color = 'blue',alpha = 0.25)
# plt.plot(master_time[three_sigma_index_err_down],master_c4_flux[three_sigma_index_err_down],
#          marker = 'x',markersize=7, color = 'red', linewidth = 0)
# plt.plot(master_time[three_sigma_index_err_up],master_c4_flux[three_sigma_index_err_up],
#          marker = 'x',markersize=7, color = 'green', linewidth = 0)
# plt.xlabel('Time (ks)', fontsize = 17)
# plt.ylabel( r'Flux  ($\dfrac{erg}{cm^2s}$) *$10^{-12}$',fontsize = 17)
# plt.title('COS G160 Observations of C IV TW Hydrae ($\Delta t = %ss$, 12.4hrs Total)' % time_resolution,fontsize = 18)


# plt.figure()
# plt.plot(master_time, master_He2_flux, color = 'blue', linewidth = 1,label = r'He II (1640)$\AA$')
# plt.fill_between(master_time, master_He2_flux-master_He2_error,master_He2_flux+ master_He2_error,color = 'blue',alpha = 0.25)
# plt.plot(master_time[three_sigma_index_err_down],master_He2_flux[three_sigma_index_err_down],
#          marker = 'x',markersize=7, color = 'red', linewidth = 0)
# plt.plot(master_time[three_sigma_index_err_up],master_He2_flux[three_sigma_index_err_up],
#          marker = 'x',markersize=7, color = 'green', linewidth = 0)
# plt.xlabel('Time (ks)', fontsize = 17)
# plt.ylabel( r'Flux  ($\dfrac{erg}{cm^2s}$) *$10^{-12}$',fontsize = 17)
# plt.title('COS G160 Observations of He II TW Hydrae ($\Delta t = %ss$, 12.4hrs Total)' % time_resolution,fontsize = 18)
