# -*- coding: utf-8 -*-
"""
This code will run inttag and calstis on each E140M observation for a desired
time resolution, as well as the 'total' time resolution

Make sure to run the following line of code before running AD_Leo_Master in
bash, in order for calstis to work

export CRDS_PATH="$HOME/crds_cache"
export CRDS_SERVER_URL="https://hst-crds.stsci.edu"
export oref="${CRDS_PATH}/references/hst/oref/"
crds bestrefs --update-bestrefs --sync-references=1 --files *.fits

"""



import os
from astropy.io import fits
import time
import stistools

#function to find all tag files in directory
def get_file_names_with_strings(str_list):
    # full_list = os.listdir(r'C:\Windows\System32\astroconda\AD_Leo_E140M')
    full_list = os.listdir('AD_Leo_E140M')
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]

    return final_list

print()
print()
print()
print('!!!Make sure you ran the CRDS procedures explained in the comments section!!!')
time.sleep(4)
print()
print()
print()
#finding all the tag files
filename_list = get_file_names_with_strings(['tag.fits'])
#setting the desired time resolution of the analysis
print()
print()
print()
print()
print('There are  ' + str(len(filename_list)) + '  files in the file list')
print('Filename List = ')
print(filename_list)
print()
print()
print()
print()
time.sleep(5)
time_resolution = int(input('Enter an integer for the time resolution of the light curves:    '))
High_Resolution = False
#now to run the inttag procedure for each file in this list
for i in filename_list:
    print('now working with  ' + i)
    print()
    print()
    time.sleep(3)
    HST_tag_file = i
    
    hdul = fits.open(HST_tag_file)

    
    TEXPTIME = hdul[0].header['TEXPTIME'] #total exposure time in seconds
    print('Total Exposure time =  %ss' % TEXPTIME)
    TEXPTIME = int(TEXPTIME)
    
    number_of_frames = int(TEXPTIME/time_resolution)
    print()
    print()
        
    filename = hdul[0].header['FILENAME'] # filename keyword
    filename_raw = filename[0:10] + 'raw_PH_' + str(time_resolution) + 's.fits'
    filename_flt = filename[0:10] + 'raw_PH_' + str(time_resolution) + 's_flt.fits'
    filename_x1d = filename[0:10] + 'raw_PH_' + str(time_resolution) + 's_x1d.fits' 
    
    filename_flt_tot = filename[0:10] + 'raw_PH_tot_' + str(TEXPTIME) + 's_flt.fits'   
    filename_raw_tot = filename[0:10] + 'raw_PH_tot_' + str(TEXPTIME) + 's.fits'
    filename_x1d_tot = filename[0:10] + 'raw_PH_tot_' + str(TEXPTIME) + 's_x1d.fits'
        
    check = os.path.exists(filename_x1d_tot)
    print('total file already exists =', check)
    if check == False:
        stistools.inttag.inttag(HST_tag_file, filename_raw_tot, rcount =1, increment = TEXPTIME,
                            highres = High_Resolution )    
            
    
    check0 = os.path.exists(filename_raw)
    print('raw file already exists =', check0)
    if check0 == True:
        print()
        # userinput = input("Overwrite " + filename_raw + '?, enter Yes or No:     ')
        userinput = 'Yes'
        if userinput == 'Yes':
            os.unlink(filename_raw)
        if userinput == 'No':
            hdul.close()
            exit()
    
    
    stistools.inttag.inttag(HST_tag_file, filename_raw, rcount =number_of_frames, increment = time_resolution,
                        highres = High_Resolution )
    hdul.close()
    stistools.calstis.calstis(filename_raw, verbose=False)
    stistools.calstis.calstis(filename_raw_tot, verbose=False)
    os.unlink(filename_flt)
    os.unlink(filename_raw)
    os.unlink(filename_flt_tot)     
    os.unlink(filename_raw_tot)
    
        
#Now I want to move all x1d files to the subfolder AD_Leo_E140M to be manually
#handled at a later time. Use the moving_files code outside of bash


    
        
        
    
    
    
    
    
    
    