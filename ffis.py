"""
Module to contain functions related to accessing FFI data
"""
import sys
from glob import glob

import numpy as np
from astropy.io import fits

############# If you move these programs you will need to update these directories and names #############
sys.path.append('/content/gdrive/My Drive/')
from . import myDir as myDir



####################################################################################
# Functions to locate a star's TESS data

def find_ffi(star, cluster):
#    homedir = os.path.expanduser("~")
    dir_ffi = myDir.project_dir(cluster)+'FFI/'
    RA = str(star["RA_ICRS"])[:7]
    DE = str(star["DE_ICRS"])[:7]
    file_ffi = glob(dir_ffi+"*"+RA+"*"+DE+"*/*.fits")
    if len(file_ffi) == 0:
        file_ffi = glob(dir_ffi+"*"+RA+"*"+DE+"*.fits")
    return(file_ffi)

def find_ffi_coord(ra, dec, cluster):
    dir_ffi = myDir.project_dir(cluster)+'FFI/'
    RA = str(ra)[:7]
    DE = str(dec)[:7]
    file_ffi = glob(dir_ffi+"*"+RA+"*"+DE+"*/*.fits")
    if len(file_ffi) == 0:
        file_ffi = glob(dir_ffi+"*"+RA+"*"+DE+"*.fits")
    return(file_ffi)



def load_ffi_fromfile(file):
    ffi_data = fits.open(file)
    images = ffi_data[1].data
    image = images[100]['Flux']
    if np.max(image) == 0:
      image = images[500]['Flux']
    return image

def load_ffi(star):
    ffi_data = fits.open(star['file_ffi'][0])
    images = ffi_data[1].data
    image = images[100]['Flux']
    if np.max(image) == 0:
      image = images[500]['Flux']
    return image

def load_ffis(star):
    ffi_data = fits.open(star['file_ffi'][0])
    images = ffi_data[1].data
    image = images['Flux']
    return image


####################################################################################
def ffi_test(ffi):
    shape = np.shape(ffi)
    val = ffi[int(shape[0]/2),int(shape[1]/2)]
    if np.isfinite(val):
      good_or_bad = 1
    else:
      good_or_bad = 0
    return good_or_bad
