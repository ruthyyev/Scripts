import montage_wrapper as montage
import numpy as np
import progressbar
import reproject
from reproject import reproject_exact
from astropy.io import fits
import os
from matplotlib import pyplot as plt
from astropy.nddata.utils import Cutout2D

os.chdir('/Users/c1541417/Documents/Edge-On/Data/SPIRE')
bar = progressbar.ProgressBar()


#Import the required table:
sources = np.loadtxt('/Users/c1541417/Documents/Edge-On/final_sample_new.csv', dtype=np.str, delimiter=',')


#Select the required columns:
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

#get the distances to each source
#NB the recessional velocities do not accurately represent distance at such
#proximity. It would be best to find their actual distances from a database
#and read them in a seperate table

for num in range(len(source_names)):
    dist = (v_list / 47.37)

#here we are finding the lin_dist which is the number of kpc per pixel for each
#map. We want to scale everything to the largest number of kpc per pixel (i.e.
#the worst resolution):

for num in bar(range(len(source_names))):

    name = source_names[num]

    for wavelength in [500]:

        header_500 = fits.open('/Users/c1541417/Documents/DustPedia/Paper/new_sample/Data/SPIRE/'+name+'_SPIRE_'+str(wavelength)+'_rotated_new.fits')
        pix_size_deg = header_500[0].header['CDELT2']
        pix_size = pix_size_deg * 3600.0
        lin_dist = ((dist_list * pix_size) / 206265.0 ) * 1000.0

#here we are finding the scaling factor to change pixel sizes to:
max_lin_dist = lin_dist[45]

print max_lin_dist
scale_factor = max_lin_dist / lin_dist
new_scale = pix_size * scale_factor


for num in range(len(source_names)):
    name = source_names[num]
    angle = pos_angle[num]
    ra = ra_list[num]
    dec = dec_list[num]
    new_pix_size = new_scale[num]

    for wavelength in [500]:

        #create a new header with the re-scaled pixel sizes:

        montage.commands.mHdr(str(ra)+', '+str(dec), width, '/Users/c1541417/Documents/Edge-On/Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.txt', pix_size = new_pix_size, rotation = angle)


#here we are applying the new headers to the fits files to create new, re-scaled
#fits files that can then be stacked:
for num in range(len(source_names)):
    name = source_names[num]
    for wavelength in [500]:
            montage.wrappers.reproject(''+name+'_SPIRE_'+str(wavelength)+'.fits', '/Users/c1541417/Documents/Edge-On/Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.fits', header='/Users/c1541417/Documents/Edge-On/Rescaled/'+name+'_SPIRE_'+str(wavelength)+'_rescaled_kpc.txt', silent_cleanup = True)
