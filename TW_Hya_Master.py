# -*- coding: utf-8 -*-
"""
costools tools webpage https://costools.readthedocs.io/en/latest/splittag.html

note: best references are in best refs file
insure all cortag and associated total x1ds are in astroconda_TW_Hya/COS_data/

Any missing best ref files can be obtained here:
    https://hst-crds.stsci.edu/

need to change 'out20' to soft coded out%s % time_resolution, and have the code
make the directory

to run, open bash and type:
    echo -e "\e[32mGreen Text\e[0m"
    conda activate cos
    cd astroconda/astroconda_TW_Hya/COS_data/
    export lref='best_refs/'
    python TW_Hya_Master_testing.py 
    
"""
import os
import shutil
import calcos
import costools
import numpy as np
import time
from astropy.io import fits


def get_file_names_with_strings(str_list):
    # full_list = os.listdir(r'C:\Windows\System32\astroconda\AD_Leo_E140M')
    full_list = os.listdir()
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]

    return final_list

def get_file_names_with_strings_PH(str_list,directory_path):
    # full_list = os.listdir(r'C:\Windows\System32\astroconda\AD_Leo_E140M')
    full_list = os.listdir(directory_path)
    final_list = [nm for ps in str_list for nm in full_list if ps in nm]

    return final_list


## moving files from downloaded folder to cwd
# to_delete_list = get_file_names_with_strings_PH(['l'],"HST/")
# print(to_delete_list)
# for i in to_delete_list:
#     to_move_list = get_file_names_with_strings_PH(['l'], "HST/%s" % i)
#     print(to_move_list)
#     for j in to_move_list:       
#         dest1 = os.getcwd()
#         source = "HST/%s/" % i
#         f = j
#         shutil.move(source+'/'+f, dest1)
            



# #deleting useless pngs
# to_delete_list = get_file_names_with_strings(['x1d.png'])
# for i in to_delete_list: os.remove(i)
counter = 1
time_resolution = 20
# filename = r"C:\Users\pahi9557\Desktop\LAB\HST_Research_Project\comps2\TWHya\COS_data\lbl208xaq_corrtag_a.fits"
# filename = 'lbl208xaq_corrtag_a.fits'
filename_list = get_file_names_with_strings(['corrtag_a'])
for i in filename_list:
    filename = i
    detector_seg = filename[-6]
    HST_ID = filename[0:9]
    print("Now working on %s" % HST_ID)
    print(counter/len(filename_list))
    print('');print('');print('');print('');print('');
    time.sleep(3)
    counter += 1
    hdul = fits.open(filename)
    # hdul.info()     #this gives an overall view of what is in each FITS file extension
    # hdul[0].header #this line prints all the header key words / info from the 1st extension
    # data = hdul[1].data
    TEXPSTRT = hdul[1].header['EXPSTART']
    TEXPEND = hdul[1].header['EXPEND'] #end time (MJD) of last exposure in the file
    TEXPTIME = hdul[1].header['EXPTIME%s' % detector_seg] #total exposure time in seconds
    TEXPTIME = int(TEXPTIME)
    print('MJD =  %ss' % TEXPSTRT)
    print('Total Exposure time =  %ss' % TEXPTIME)
    Good_time_interval = hdul[2].data; Good_time_interval = np.asarray(Good_time_interval[0]);
    hdul.close()
    
    print()
    print()
    # time_resolution = int(input('Enter an integer for the time resolution of the light curves:    '))
    
    number_of_frames = int(TEXPTIME/time_resolution)
    print()
    print()
    print('Number of frames = %s' % number_of_frames)
    time.sleep(3)
    #split_tagging the a segment
    print('Now completing %s' % filename)
    costools.splittag.splittag(filename, HST_ID + '_split_%s' % detector_seg,
            starttime=0., increment=time_resolution,endtime = number_of_frames*time_resolution,
            time_list="")
    #split_tagging the b segment
    filename = list(filename)
    filename[-6] = 'b'
    filename = "".join(filename)
    detector_seg = filename[-6]
    print('Now completing %s' % filename)
    costools.splittag.splittag(filename, HST_ID + '_split_%s' % detector_seg,
            starttime=0., increment=time_resolution,endtime = number_of_frames*time_resolution,
            time_list="")
    
    os.makedirs("out20/%s/" % HST_ID)
    # #make x1ds
    to_calcos_list = get_file_names_with_strings(['split'])
    for i in to_calcos_list:
        calcos.calcos(i, outdir="out20/%s/" % HST_ID)
        #deleting needless files from sub-directory
        to_delete_list = get_file_names_with_strings_PH(['tra','tag', 'counts','flt'],"out20/%s/" % HST_ID)
        for i in to_delete_list: os.remove("out20/%s/%s" % (HST_ID,i))
    #deleting needless files from main directory
    to_delete_list = get_file_names_with_strings(['split'])
    for i in to_delete_list: os.remove(i)
    x1d_list = get_file_names_with_strings(['_x1d'])
    for i in x1d_list:
        new_name = i[0:9] + '-total-x1d.fits'
        os.rename(i,new_name)
    #move total x1d file into sub direcotry
    source = os.getcwd()
    dest1 = "out20/%s/" % HST_ID
    f = HST_ID + '-total-x1d.fits'
    shutil.copy(source+'/'+f, dest1)

print('completed %s number of files' % counter)