from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt

#Import the required table:
sources = np.loadtxt('/Users/c1541417/Documents/DustPedia/Paper/new_sample/Sample/final_sample.csv', dtype=np.str, delimiter=',')


#Select the required columns, i.e. source names and position angles:
source_names = sources[1:,0]


for name in source_names:
    for wavelength in [500]:
        try:
            image_list = [fits.getdata('/Users/c1541417/Documents/DustPedia/Paper/new_sample/kpc_scaling/Cutouts/'+name+'_'+str(wavelength)+'_cutout.fits')]
        except:
            print 'could not find '+name+''

image_concat = []
for image in image_list:
    image_concat.append(image)


final_image = np.sum(image_concat, axis=0)

plt.imshow(final_image)


outfile = '/Users/c1541417/Documents/DustPedia/Paper/new_sample/kpc_scaling/Stacked/500_cutout_stacked.fits'

hdu = fits.PrimaryHDU(final_image)
hdu.writeto(outfile, clobber=True)
