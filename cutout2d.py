from astropy.nddata.utils import Cutout2D
import numpy as np
from astropy.io import fits
import os
from matplotlib import pyplot as plt

#set the working directory
os.chdir('/Users/c1541417/Documents/DustPedia/Paper/new_sample/kpc_scaling/Rescaled/')

#Import the required table:
sources = np.loadtxt('/Users/c1541417/Documents/DustPedia/Paper/new_sample/Sample/final_sample.csv', dtype=np.str, delimiter=',')

#Select the required columns, i.e. source names and position angles:
source_names = sources[1:,0]

#make empty lists to append information to
x_position = []
y_position = []
cutouts = []

#set a size for the cutouts in pixels:
size = 30

for num in range(len(source_names)):
    name = source_names[num]
    data_500 = fits.open(''+name+'_SPIRE_500_rescaled_kpc.fits')

    #get the galaxy central pixels from the header files
    x_coords = data_500[0].header['CRPIX1']
    y_coords = data_500[0].header['CRPIX2']

    #append information to empty lists
    x_position.append(x_coords)
    y_position.append(y_coords)

    x_pos = x_position[num]
    y_pos = y_position[num]

    data = data_500[0].data

    #make the cutouts
    cutouts = Cutout2D(data, (x_pos, y_pos), (10, 25))

    #save the resulting cutouts as fits files
    fits.writeto('/Users/c1541417/Documents/DustPedia/Paper/new_sample/kpc_scaling/Cutouts/'+name+'_500_cutout.fits', cutouts.data)
