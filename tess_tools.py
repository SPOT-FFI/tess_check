# Plotting program
import matplotlib.pyplot as plt

# Search for files
from glob import glob

# Good for reading/writing data tables
import pandas as pd

# Better math, numbers, and array functions
import numpy as np

from astropy.timeseries import LombScargle
from astropy.io import fits


# You have to give permission to the notebook to access files stored in your Google Drive.
from google.colab import drive
drive.mount('/content/gdrive',force_remount=True)

import os
import sys
############# If you move these programs you will need to update these directories and names #############
sys.path.append('/content/gdrive/Shareddrives/')
from tess_check import myDir as myDir


def findel(val, array):
	array = np.array(array)
	adiff = np.abs(array-val)
	id = np.argmin(adiff)
	return id

# Download TESS
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
import astropy.units as u


# CPM modules
from google.colab import drive
drive.mount('/content/gdrive',force_remount=True)
sys.path.append('/content/gdrive/Shareddrives/DouglasGroup/tess_check/')
from tess_check import tess_cpm as tess_cpm
from tess_check import tesscheck as tesscheck

# SAP modules
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.stats import sigma_clipped_stats

# Make sheet
def make_sheet(project, file):
    # read file
    # organize data
    # generate sheet

    columns = ['DR2Name','RA','Dec','Gmag','BP_RP','Prot','Prot_LS','Power_LS','TESS_Data','Notes']



# Load sheet
def load_sheet(project_name):
    from google.colab import auth
    auth.authenticate_user()
    import gspread
    from oauth2client.client import GoogleCredentials
    gc = gspread.authorize(GoogleCredentials.get_application_default())

    dir_project = myDir.project_dir(project_name)
    sheet_id_file = os.path.join(dir_project,f"{project_name}.txt")

    f = open(sheet_id_file,"r")
    doc_id = f.readline()
    f.close()
    sheet = gc.open_by_key(doc_id)
    return sheet

def download_tess(sample_name,mag_cut=16.0,gbr_cut=0.5):
    # Check for directories
    # Main level
    dir_project = myDir.project_dir(sample_name)
    levels = ['','/FFI','/CPM','/Plots','/Panels','/SAP']

    for level in levels:
        path_level = dir_project+level
        dir_level = glob(path_level)
        if np.size(dir_level) == 0:
            os.mkdir(path_level)
            print('Created Dir: '+path_level)

    dir_ffi = dir_project + 'FFI/'

    # Load list
    worksheet = load_sheet(sample_name)
    target_table = worksheet.worksheet('Targets')
    target_data_val = target_table.get_all_values()
    target_data = pd.DataFrame.from_records(target_data_val[1:],columns=target_data_val[0])

    gbr = target_data['BP_RP'].to_numpy()
    gmag = target_data['Gmag'].to_numpy()

    # Trim list

    no_gbr = np.where(gbr == 'nan') # 20 stars with no gbr?
    # if there isn't a color, then download anyway if it passes the mag_cut

    number_of_stars = len(gbr)
    print(number_of_stars)
    for i in range(number_of_stars):
        if (gbr[i] != 'nan') & (gmag[i] != 'nan') & (len(target_data['RA'][i]) != 0):
            if (float(gbr[i])>gbr_cut) & (float(gmag[i])<mag_cut):
                # test if files exist
                RA_test = str(np.round(float(target_data['RA'][i]),6))
                DE_test = str(np.round(float(target_data['Dec'][i]),6))
                file_ffi = glob(dir_ffi+"*"+RA_test+"*"+DE_test+"*/*.fits")
                if len(file_ffi) == 0:
                    file_ffi = glob(dir_ffi+"*"+RA_test+"*"+DE_test+"*.fits")
                if np.size(file_ffi) == 0:
                    cutout_coord = SkyCoord(target_data['RA'][i], target_data['Dec'][i],unit = "deg")
                    manifest = Tesscut.download_cutouts(cutout_coord, size=40, path=dir_ffi)
                    print(str(i+1)+' of '+str(number_of_stars)+' stars')
    print('Done')


    # Check data sizes
    # delete any small files
    # Stephanie has deleted this because all of the NGC 2632 files were <100b
#    files = glob(dir_ffi+'*')
#    for file in files:
#        bsize = os.path.getsize(file) #size in bytes
#        if (bsize < 200):
#            os.remove(file)

# Define tools

def cpm_doplot(file,t_row=20, t_col=20, excl=4):
    cpm = tess_cpm.CPM(file, remove_bad=True)
    cpm.set_target(t_row, t_col)
    cpm.set_exclusion(excl)
    cpm.set_predictor_pixels(256, method='similar_brightness') # cosine_similarity
    cpm.lsq(0.1, rescale=True, polynomials=False)
    tess_cpm.summary_plot(cpm, 15)
    #aperture_lc, lc_matrix = cpm.get_aperture_lc(box=1, show_pixel_lc=True, show_aperture_lc=True)
    time = cpm.time
    flux = cpm.rescaled_target_fluxes - cpm.lsq_prediction
    return time, flux

def cpm_do(file,t_row=20, t_col=20, excl=4):
    cpm = tess_cpm.CPM(file, remove_bad=True)
    cpm.set_target(t_row, t_col)
    cpm.set_exclusion(excl)
    cpm.set_predictor_pixels(256, method='similar_brightness') # cosine_similarity
    cpm.lsq(0.1, rescale=True, polynomials=False)
#    tess_cpm.summary_plot(cpm, 15)
    #aperture_lc, lc_matrix = cpm.get_aperture_lc(box=1, show_pixel_lc=True, show_aperture_lc=True)
    time = cpm.time
    flux = cpm.rescaled_target_fluxes - cpm.lsq_prediction
    return time, flux

def cpm_sector(file):
    start = file.find('tess-s')
    sector = file[start+8:start+8+2]
    return sector

def cpm_run(files,t_row=20, t_col=20, excl=4):
    time, flux, sector = [], [], []
    for f in files:
      # test FFI image
      ffi_image = tesscheck.load_ffi_fromfile(f)
      if tesscheck.ffi_test(ffi_image) == 1:
        this_time, this_flux = cpm_do(f,t_row=t_row,t_col=t_col,excl=excl)
        time.append(this_time)
        flux.append(this_flux)
        this_sector = cpm_sector(f)
        sector_list = np.zeros(len(this_time)) + np.int(this_sector)
        sector.append(sector_list)

    if np.size(flux) == 0:
      print('light curve not available')
      return -1,-1,-1
    time = np.concatenate([t for t in time])
    flux = np.concatenate([b for b in flux])
    sector = np.concatenate([s for s in sector])
    # double check light curve for NANs
    ifin = np.where(np.isfinite(flux) == True)
    if np.size(ifin) == 0:
      print('light curve not available')
      return -1,-1,-1
    else:
      time = time[ifin[0]]
      flux = flux[ifin[0]]
      sector = sector[ifin[0]]

    return time, flux, sector

def cpm_per(time, flux):
    fmax = 1/0.1 # 0.1 days
    fmin = 1/50 # 50 days
    freq, power = LombScargle(time, flux).autopower(minimum_frequency=fmin, maximum_frequency=fmax)
    return 1/freq[np.argmax(power)],np.max(power)

def cpm_perplot(time, flux):
    fmax = 1/0.1 # 0.1 days
    fmin = 1/50 # 50 days
    freq, power = LombScargle(time, flux).autopower(minimum_frequency=fmin, maximum_frequency=fmax)

    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(14,6))
#    ax1.scatter(time, flux, markersize=4)
    ax1.scatter(time, flux,s=3)
    ax1.set_ylim(-0.05,0.05)
    ax2.plot(1/freq, power)
    ax1.set_title('CPM Light Curve', fontsize=20)
    ax2.set_title('Lomb-Scargle', fontsize=20)
    plt.show()
    per = 1/freq[np.argmax(power)]
    print(per)
    return per

def make_cpm(sample_name):
    # Main level
    dir_project = myDir.project_dir(sample_name)
    dir_ffi = dir_project + 'FFI/'
    dir_cpm = dir_project + 'CPM/'

    # Load list
    sheet_name =  sample_name
    from google.colab import auth
    auth.authenticate_user()
    import gspread
    from oauth2client.client import GoogleCredentials
    gc = gspread.authorize(GoogleCredentials.get_application_default())
    # load the table
    worksheet = load_sheet(sample_name)

    target_table = worksheet.worksheet('Targets')
    target_data_val = target_table.get_all_values()
    target_data = pd.DataFrame.from_records(target_data_val[1:],columns=target_data_val[0])

    gbr = target_data['BP_RP'].to_numpy()
    gmag = target_data['Gmag'].to_numpy()

    number_of_stars = len(gbr)

    for i in range(number_of_stars):
        csv_file = dir_cpm+'GaiaDR2_'+target_data['DR2Name'][i]+'.csv'
        lc_test = glob(csv_file)
        if np.size(lc_test) == 0:
            RA_test = str(target_data['RA'][i])[:7]
            DE_test = str(target_data['Dec'][i])[:7]
            file_ffi = glob(dir_ffi+"*"+RA_test+"*"+DE_test+"*/*.fits")
            if len(file_ffi) == 0:
                file_ffi = glob(dir_ffi+"*"+RA_test+"*"+DE_test+"*.fits")
            if np.size(file_ffi) != 0:
                time, flux, sector = cpm_run(file_ffi,t_row=20, t_col=20, excl=5)
                if np.size(time)>1:
                    df = pd.DataFrame({'time' : time, 'flux' : flux, 'sector' : sector})
                    df.to_csv(csv_file, index=False)
                    print('Star '+str(i+1)+' of '+str(number_of_stars))


def load_ffis(star):
    ffi_data = fits.open(star['file_ffi'][0])
    images = ffi_data[1].data
    times = images['Time']
    image = images['Flux']
    quality = images['Quality']
    iq = np.where(quality == 0)
    good_times = times[iq[0]]
    good_images = image[iq[0]]

    return good_times, good_images

def tess_sector(file):
    start = file.find('tess-s')
    sector = file[start+8:start+8+2]
    return sector

def make_sap(project_name):

    # Main level
    dir_project = myDir.project_dir(project_name)
    dir_ffi = dir_project + 'FFI/'
    dir_sap = dir_project + 'SAP/'

    # We need to give Colab access to the data stored in our Drive.
    from google.colab import drive
    drive.mount('/content/gdrive/', force_remount=True)

    # Load list
    worksheet = load_sheet(project_name)
    target_table = worksheet.worksheet('Targets')
    target_data_val = target_table.get_all_values()
    target_data = pd.DataFrame.from_records(target_data_val[1:],columns=target_data_val[0])

    RA = target_data['RA'].to_numpy()
    Dec = target_data['Dec'].to_numpy()

    number_of_stars = len(target_data)
    print('Number of stars: '+str(number_of_stars))

    for i in range(number_of_stars):
        csv_file = dir_sap+'GaiaDR2_'+target_data['DR2Name'][i]+'-SAP.csv'
        lc_test = glob(csv_file)
        if (np.size(lc_test) == 0) & (RA[i] != 'nan'):
            RA_test = str(RA[i])[:7]
            DE_test = str(Dec[i])[:7]
            file_ffi = glob(dir_ffi+"*"+RA_test+"*"+DE_test+"*/*.fits")
            if len(file_ffi) == 0:
                file_ffi = glob(dir_ffi+"*"+RA_test+"*"+DE_test+"*.fits")
            if np.size(file_ffi) != 0:
                time, flux, sector = sap_run(file_ffi,t_row=20, t_col=20, rad=4, sky=(6,15))
                if np.size(time)>1:
                    df = pd.DataFrame({'time' : time, 'flux' : flux, 'sector' : sector})
                    df.to_csv(csv_file, index=False)
                    print('Star '+str(i+1)+' of '+str(number_of_stars))



def sap_do(file, t_row=20, t_col=20, rad=4, sky=(6,15)):
  # load the fits file
  ffi_data = fits.open(file)
  images = ffi_data[1].data
  times = images['Time']
  image = images['Flux']
  quality = images['Quality']
  iq = np.where(quality == 0)
  good_times = times[iq[0]]
  good_images = image[iq[0]]
  # prep for photometry
  number_of_obs = len(good_times)
  mags = np.zeros(number_of_obs)
  # loop over observations
  for i_image in range(number_of_obs):
    ffi_image = good_images[i_image]
    data = ffi_image-np.min(ffi_image)
    positions = [(t_row,t_col)]
    aperture = CircularAperture(positions, r=rad)
    annulus_aperture = CircularAnnulus(positions, r_in=sky[0], r_out=sky[1])
    annulus_masks = annulus_aperture.to_mask(method='center')
    annulus_data = annulus_masks[0].multiply(data)
    annulus_data_1d = annulus_data[annulus_masks[0].data > 0]
    _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
    bkg_median = np.array(median_sigclip)
    phot = aperture_photometry(data, aperture)
    phot['annulus_median'] = bkg_median
    phot['aper_bkg'] = bkg_median * aperture.area
    phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
    mags[i_image] = float(phot['aper_sum_bkgsub'])

  mags /= np.median(mags)
  return good_times, mags

def sap_run(files, t_row=20, t_col=20, rad=4, sky=(6,15)):
    time, flux, sector = [], [], []
    for f in files:
      # test FFI image
      ffi_image = tesscheck.load_ffi_fromfile(f)
      if tesscheck.ffi_test(ffi_image) == 1:
        this_time, this_flux = sap_do(f, t_row=t_row, t_col=t_col, rad=4, sky=(6,15))
        time.append(this_time)
        flux.append(this_flux)
        this_sector = tess_sector(f)
        sector_list = np.zeros(len(this_time)) + np.int(this_sector)
        sector.append(sector_list)

    if np.size(flux) == 0:
      print('light curve not available')
      return -1,-1,-1
    time = np.concatenate([t for t in time])
    flux = np.concatenate([b for b in flux])
    sector = np.concatenate([s for s in sector])
    # double check light curve for NANs
    ifin = np.where(np.isfinite(flux) == True)
    if np.size(ifin) == 0:
      print('light curve not available')
      return -1,-1,-1
    else:
      time = time[ifin[0]]
      flux = flux[ifin[0]]
      sector = sector[ifin[0]]
      flux /= np.median(flux)
    return time, flux, sector


def list_to_string(list,delim=' '):
  list_as_string = list[0]
  for i in range(len(list)-1):
    id = i+1
    list_as_string += (delim + list[id])
  return list_as_string
