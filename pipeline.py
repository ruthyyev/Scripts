#PIPELINE.PY
#A script to read in data files, rotate and regrid them to the same no. of kpc
#per pixel. It then makes cutouts of each galaxy based on their WCS info and
#stacks the cutouts to produce one map with all stacked galaxies. -- NB THE
#FINAL IMAGE IS NOT ACTUALLY STACKED!! - sort your shit out ruth

#-------------------------------------------------------------------------------
#IMPORT REQUIRED MODULES:
#-------------------------------------------------------------------------------
import montage_wrapper as montage
import numpy as np
import progressbar
import reproject
from reproject import reproject_exact
from astropy.io import fits
import os
from matplotlib import pyplot as plt
from astropy.nddata.utils import Cutout2D

#Set the current working directory:
os.chdir('/Users/c1541417/Documents/Edge-On/')
#Include progress bar when running the script
bar = progressbar.ProgressBar()


#Import the required table and select the required columns:
sources = np.loadtxt('final_sample_new.csv', dtype=np.str, delimiter=',')
#source names:
source_names = sources[1:,0]
#position angles:
pos_angle0 = sources[1:,65]
#the right ascentions of the sources:
ra_list = sources[1:,1]
#the declinations of the sources:
dec_list = sources[1:,2]
#the recessional velocities of the sources:
v_0 = sources[1:,5]
#the distance to the targets
dist_0 = sources[1:,66]

#change the velocities to floats as calculations won't work with strings!
v_list = v_0.astype(np.float)
dist_list = dist_0.astype(float)

#adjust the position angles as such so that they will all rotate to the same
#angle:
pos_angle1 = pos_angle0.astype(np.float)
pos_angle = [360. - x for x in pos_angle1]

#this is the width of the image that montage outputs:
width = 5.

print source_names[43]


#-------------------------------------------------------------------------------
#RE-SCALE THE IMAGES TO THE SAME KPC PER PIXEL
#-------------------------------------------------------------------------------
#here we are finding the lin_dist which is the number of kpc per pixel for each
#map. We want to scale everything to the largest number of kpc per pixel (i.e.
#the worst resolution):
wave = [250, 350, 500]
"""
for num in range(len(source_names)):
    dist = (v_list / 47.37)
    name = source_names[num]

#SORT THIS SHIT OUT!:
    header_250 = fits.open('Data/SPIRE/'+name+'_SPIRE_250.fits')
    pix_size_deg_250 = header_250[0].header['CDELT2']
    pix_size_250 = pix_size_deg_250 * 3600.0
    lin_dist_250 = ((dist_list * pix_size_250) / 206265.0 ) * 1000.0

    header_350 = fits.open('Data/SPIRE/'+name+'_SPIRE_350.fits')
    pix_size_deg_350 = header_350[0].header['CDELT2']
    pix_size_350 = pix_size_deg_350 * 3600.0
    lin_dist_350 = ((dist_list * pix_size_350) / 206265.0 ) * 1000.0

    header_500 = fits.open('Data/SPIRE/'+name+'_SPIRE_500.fits')
    pix_size_deg_500 = header_500[0].header['CDELT2']
    pix_size_500 = pix_size_deg_500 * 3600.0
    lin_dist_500 = ((dist_list * pix_size_500) / 206265.0 ) * 1000.0

#here we are finding the scaling factor to change pixel sizes to:


max_lin_dist = lin_dist[43]
scale_factor = max_lin_dist / lin_dist
new_scale = pix_size * scale_factor
print new_scale[43]


for num in range(len(source_names)):
    name = source_names[num]
    angle = pos_angle[num]
    ra = ra_list[num]
    dec = dec_list[num]
    new_pix_size = new_scale[num]

    for wavelength in [250, 350, 500]:

        #create a new header with the re-scaled pixel sizes:

        montage.commands.mHdr(str(ra)+', '+str(dec), width,
                              'Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.txt',
                              pix_size = new_pix_size, rotation = angle)


#here we are applying the new headers to the fits files to create new, re-scaled
#fits files that can then be stacked:
for num in range(len(source_names)):
    name = source_names[num]
    for wavelength in [250, 350, 500]:
            montage.wrappers.reproject(''+name+'_SPIRE_'+str(wavelength)+'.fits',
                                       '/Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.fits',
                                       header='Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.txt',
                                       silent_cleanup = True)
"""
#-------------------------------------------------------------------------------
#MAKE CUT-OUTS OF THE DATA - THE GALAXY AT THE CENTRE
#-------------------------------------------------------------------------------

#make empty lists to append information to
x_position = []
y_position = []
cutouts = []

for num in range(len(source_names)):
    for wavelength in [250]:
        name = source_names[num]
        data = fits.open('Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.fits')

    #get the galaxy central pixels from the header files
    x_coords = data[0].header['CRPIX1']
    y_coords = data[0].header['CRPIX2']

    #append information to empty lists
    x_position.append(x_coords)
    y_position.append(y_coords)

    x_pos = x_position[num]
    y_pos = y_position[num]

    fin_data = data[0].data

    #make the cutouts
    cutouts = Cutout2D(fin_data, (x_pos, y_pos), (50, 50), mode = 'partial', fill_value = 0.)

    #save the resulting cutouts as fits files
    fits.writeto('Cutouts/'+name+'_'+str(wavelength)+'_cutout.fits',
                     cutouts.data, clobber = True)


#-------------------------------------------------------------------------------
#STACK THE IMAGES
#-------------------------------------------------------------------------------

outputs = []

for wavelength in [250, 350, 500]:

    image_list = []

    for name in source_names:

        try:
            image_list.append(fits.getdata('Cutouts/'+name+'_'+str(wavelength)+'_cutout.fits'))
        except:
            print 'could not find '+name+''

    image_list = [i for i in image_list if i.shape == (50, 50)]
    image_list = [np.ma.masked_invalid(i) for i in image_list]

#    for im in image_list:
#        print(im.shape)

    try:
        outputs.append(np.ma.sum(image_list, axis=0))
    except:
        pass

    plt.figure()
    plt.imshow(outputs[-1])

    outfile = 'Results/'+str(wavelength)+'_cutout_stacked.fits'
    fits.writeto(outfile, outputs[-1].data, clobber = True)
