# this will be the main program for inspecting TESS light curves for stellar rotation

# Import relevant modules
#%matplotlib inline
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
#import matplotlib.cm as cm
#import matplotlib
import matplotlib.gridspec as gridspec
#from astropy.visualization import astropy_mpl_style
from glob import glob
from astropy.io import fits
import warnings
warnings.filterwarnings('ignore')
#from astroquery.vizier import Vizier
#from astropy.coordinates import SkyCoord
#import astropy.units as u
from astropy.timeseries import LombScargle

import os
import sys
############# If you move these programs you will need to update these directories and names #############
sys.path.append('/content/gdrive/My Drive/')
from tess_check import myDir as myDir

import time
####################################################################################

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


def list_to_string(list,delim=' '):
  list_as_string = list[0]
  for i in range(len(list)-1):
    id = i+1
    list_as_string += (delim + list[id])
  return list_as_string

# Auto run
def tesscheck_auto(project_name, tess_cycle=1, redo=False):
    #cluster = 'ComaBer'
    user = 'Auto'
    target_table = get_prot_table('Auto',project_name)
    target_data_val = target_table.get_all_values()
    target_data = pd.DataFrame.from_records(target_data_val[1:],columns=target_data_val[0])

    if tess_cycle == 1:
        cycle_sectors = ['1','2','3','4','5','6','7','8','9','10','11','12','13']
    if tess_cycle == 2:
        cycle_sectors = ['14', '15', '16','17','18','19','20','21','22','23','24','25','26']

    print('Assembling target list...')
    #star_list = stars_todo(target_table)
    star_list = stars_todo(target_data)
    number_stars = len(star_list)
    print(str(number_stars)+' stars to analyze')
    if number_stars == 0:
        print('no stars remaining.')
    else:
        for i in range(number_stars):
            if star_list[i][0] != 'n':
                star = make_star(target_data,star_list[i])
                star = make_star(target_data,star_list[i])
                print(str(i)+'  '+star_list[i])
                tstar = initiate_star(star,project_name,user=user)
                tstar['which_sectors'] = cycle_sectors
                display_tess_lite_v2(tstar, save = False, noplot = True)
                update_prot_table(target_table, tstar)
    return

####################################################################################
# Main program
def tesscheck_run_v1():
	# Identify the user
	user = tess_user()

	# Load the Prot table
	prot_table = get_prot_table(user)
	
#	rows = prot_table.get_all_values()
#	print(rows)

	cluster = 'NGC_7092'
	file = glob('/content/gdrive/My Drive/Tables/'+cluster+'-Catalog.csv')
	clu = pd.read_csv(file[0])
	gmag = clu['GMAG']
	bprp = clu['BP_RP']
	RA = clu['RA_ICRS']
	Dec = clu['DE_ICRS']
#	star = clu.iloc[121]
	star = clu.iloc[276]
	tstar = initiate_star(star,cluster,user=user)
	#display_tess_lite(tstar, save = True)
	
	display_tess_lite_v2(tstar, save = False, noplot = True)
	
	update_prot_table(prot_table, tstar, user)
	
	return tstar


####################################################################################
# Identify the User
def tess_user(project_name):
    #dir_project = project_dir(project_name)
    status = read_status(project_name)

    #file_status = glob(dir_project + 'Status.txt')
    #file_open = open(file_status[0], "r")
    #lines = file_open.readlines()
    #user_line = lines[8]
    users = status['Users'].split(' ')
    #users = users[1:-1]
    number_of_users = len(users)
    
    print('Which user? Press...')
    for i in range(number_of_users):
        print('   '+str(i+1)+' for '+users[i])
    print('   '+str(i+2)+' for Other')

# while loop

#  val = input("Enter number for your name: ") 
    val = input()
    user = ''
    #user_found = 0
    
    if val.isdigit() == False:
        print('No user selected.')
        user = 'None'
        return user
    else: 
    # is the number in the range?
        if (float(val) > (number_of_users+1)) | (float(val) == 0):
            print('Out of bounds')
            user = 'Noone'
            return user
        if (float(val) <= (number_of_users)) & (float(val) != 0):
            id = int(val)
            user = users[id-1]
            print(user + ' is logged in.')
            add_user_to_sheet(project_name, user)
            return user
        if float(val) == (number_of_users+1):
            print('Other selected. Need to make sheet and update Status.txt')
            print('Other: What name?')
            other_name = input()
            other_name = other_name.replace(" ", "")
            other_name = other_name.lower()
            other_name = other_name.capitalize()
            #make_status(project_name, add_user=other_name, add_step = '', change_auto='', reset=False):
            user = other_name # prompt for it
            iu = np.where(np.array(users) == user)
            if np.size(iu) > 0: 
                print('User already exists, logging in.')
                add_user_to_sheet(project_name, user)
                return user
            else:
                make_status(project_name, add_user=user, add_step = '', change_auto='', reset=False)
                add_user_to_sheet(project_name, user)
                print(user + ' profile is created and user is logged in.')
            return user

####################################################################################
def get_prot_table(user,project_name):
	sheet_name = project_name
	from google.colab import auth
	auth.authenticate_user()
	import gspread
	from oauth2client.client import GoogleCredentials
	gc = gspread.authorize(GoogleCredentials.get_application_default())
	# load the table
	worksheet = load_sheet(sheet_name)
	if user == 'Auto':
		table = worksheet.worksheet('Targets')
	else:
		table = worksheet.worksheet(user)
	return table
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


def find_sap(star, cluster):
    dir_sap = myDir.project_dir(cluster)+'SAP/'
    if star["DR2NAME"][0] == 'G':
        star_name = star["DR2NAME"][9:]
    else:
        star_name = star["DR2NAME"]
    file_sap = glob(dir_sap + "*" + star_name + "*.csv")
    return(file_sap)

def find_cpm(star, cluster):
    dir_cpm = myDir.project_dir(cluster)+'CPM/'
    if star["DR2NAME"][0] == 'G':
        star_name = star["DR2NAME"][9:]
    else:
        star_name = star["DR2NAME"]
    file_cpm = glob(dir_cpm + "*" + star_name + "*.csv")
    return(file_cpm)

def load_cpm_fromfile(file):
    lc = pd.read_csv(file)
    return lc

def load_cpm(star):
    lc = pd.read_csv(star['file_cpm'][0])
    return lc

def load_sap(star):
    lc = pd.read_csv(star['file_sap'][0])
    return lc

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
def make_star(target_data, dr2_now):
    star = {'DR2NAME': '',
            'RA_ICRS': 0.,
            'DE_ICRS': 0.,
            'GMAG': 0., 
            'BP_RP': 0.}

    iloc = np.where(dr2_now == target_data['DR2Name'])
    id = iloc[0][0]

    star['DR2NAME'] = dr2_now
    star['RA_ICRS'] = target_data['RA'][id]
    star['DE_ICRS'] = target_data['Dec'][id]
    star['GMAG'] = target_data['Gmag'][id]
    star['BP_RP'] = target_data['BP_RP'][id]

    return star
####################################################################################
def initiate_star(star, cluster, user='NONE', blank=False):
    star_data = {'User' : user,
                'Source': '',
                'RA':0.,'Dec':0.,
                'Gmag':0., 'gbr':0.,
                #'Cluster':star['Cluster'],
                 'Cluster':cluster,
                # data files
                'file_ffi':'', 'file_cpm':'', 'file_sap':'', 'file_cdips':'',
                'exist_ffi': 0, 'exist_cpm':0, 'exist_cdips':0, 'exist_sap':0,
                'SAP_line':False,
                'number_sectors':0,
                'which_sectors':[''],
                'sector_list':[''],
                # lc arrays
                # LC option and Period results
                'which_LC':'CPM', # default is CPM
                'Prot_LS':0., 'is_it_double':0, # the Lomb-Scargle period and if it should be doubled
                'Power_LS':0., # the Lomb-Scargle period and if it should be doubled
                'Prot_final':0.,
                'Amplitude':0.,
                'Multi':0, # if it is multi, then set this to multi=1
                'Flares':0,
                'Notes':'', # anything noteworthy about this star?
                'LC_Quality':1, # 1 = modulate, 0 is flat, -1 is garbage
                'LC_Action':'', # 
                # Plotting options
                'x_min':0.0,'x_max':0.0, # time range
                'y_min':0.0,'y_max':0.0, # flux range, will be calculated during first iteration, and then adjusted by user
                # LS options
                'pmin':0.1,'pmax':30., # default period range for LS analysis
                'pxlog':0
                }

    if blank == True:
        return star_data

    star_data['Source'] = star["DR2NAME"]
    star_data['RA'] = star["RA_ICRS"]
    star_data['Dec'] = star["DE_ICRS"]
    star_data['Gmag'] = star["GMAG"]
    star_data['gbr'] = star["BP_RP"]

    # Once the blank data dictionary is created, test/load data into it
    exist = np.array([0,0,0,0])
    file_ffi = find_ffi(star,cluster)
    file_cpm = find_cpm(star,cluster)
    file_sap = find_sap(star,cluster)
    file_cdips = ''
    #file_cdips = find_cdips(star)
    if len(file_ffi) > 0:
        exist[0] = 1
        # star_data['ffi_image'] = load_ffi(file_ffi)
        star_data['file_ffi'] = file_ffi
        star_data['exist_ffi'] = 1
    if len(file_cpm) > 0:
        exist[1] = 1
        #lc_cpm = load_cpm(file_cpm)
        star_data['file_cpm'] = file_cpm            
        star_data['exist_cpm'] = 1
    if len(file_sap) > 0:
        exist[2] = 1
        #lc_cpm = load_cpm(file_cpm)
        star_data['file_sap'] = file_sap
        star_data['exist_sap'] = 1
    if len(file_cdips) > 0:
        exist[3] = 1
        #lc_cdips = load_cdips(file_cdips)
        star_data['file_cdips'] = file_cdips
        star_data['exist_cdips'] = 1
    if exist.sum() == 0:
        print('No data for this star')
        return star_data
    else:
        return star_data
####################################################################################
# modified version of the display program. needs to be copied into tesscheck.py
def display_tess_lite_v2(tstar,save = False,noplot = False):
    time_1 = time.time()
    # plotting defaults
    axis_fontsize = 16
    import matplotlib.pylab as pylab
    params = {'axes.labelsize': 16,'axes.titlesize': 16,'xtick.labelsize': 14,'ytick.labelsize': 14}
    pylab.rcParams.update(params)
    
#cpm data for star
    if tstar['exist_cpm'] == 0:
        return
    
    if tstar['which_LC'] == 'CPM':
        lc_cpm = load_cpm(tstar)
#        if tstar['which_sectors'] != 'All':
#            ifin = np.where((np.isfinite(lc_cpm['flux']) == True) & (lc_cpm['sector'] == int(tstar['which_sectors'])))
#            if np.size(ifin) == 0:
#                print('sector not available, reverting back to All')
#                ifin = np.where(np.isfinite(lc_cpm['flux']) == True)
#                tstar['which_sectors'] = 'All'
        time_all = lc_cpm['time']
        flux_all = lc_cpm['flux']
        sector_all = lc_cpm['sector']
        lc_title = 'TESS Light Curve (Calibrated with Causal Pixel Modeling)'
    if tstar['which_LC'] == 'SAP':
        lc_sap = load_sap(tstar)
        time_all = lc_sap['time']
        flux_all = lc_sap['flux']
        flux_all /= np.nanmedian(flux_all)
        flux_all -= 1
        sector_all = lc_sap['sector']
        lc_title = 'TESS Light Curve (Extracted with Simple Aperture Photometry)'

    # what if we say just load the whle thing, whether SAP or CPM, then we handle the sectors...
    ifin = np.where(np.isfinite(flux_all)) #find infs
    unique_sectors = np.unique(sector_all).astype(int) #all the sectors for each data point from table
    unique_sectors = list(map(str,unique_sectors)) #unique instances
    # save into 
    tstar['sector_list'] = unique_sectors #unique sectors saved as a list
    length_all = len(flux_all)
    use_these = np.zeros(length_all)
    length_which_sectors = len(tstar['which_sectors']) #tstars['which_sectors'] is blank until this runs once,
    if tstar['which_sectors'] != ['']: #skips this on first run when it's blank
        for index_sectors in range(length_which_sectors):
          id_sector_match = np.where((float(tstar['which_sectors'][index_sectors]) == sector_all) & (np.isfinite(flux_all) == True))
          if len(id_sector_match[0]) > 0:
            use_these[id_sector_match[0]] = 1

    ifin = np.where(use_these == 1)
    if len(ifin[0]) == 0:
      ifin = np.where(np.isfinite(flux_all) == True)
      print('all points ruled out, reverting to all points')
      use_these = np.zeros(length_all)
      use_these[ifin[0]] = 1

    if float(tstar['y_max'])>0:
#        print('trimming')
        #ifin = np.where((flux>float(tstar['y_min'])) & (flux<float(tstar['y_max'])))
        iyra = np.where(abs(100*flux_all)>float(tstar['y_max']))
        if len(iyra[0]) > 0:
          use_these[iyra[0]] = 0

    if ((tstar['x_min'] != 0) | (tstar['x_max'] != 0)):
        ixra = np.where((time_all-min(time_all) < float(tstar['x_min'])) | (time_all-min(time_all) > float(tstar['x_max'])))
        if len(ixra[0]) > 0:
          use_these[ixra[0]] = 0

 #    ifin = np.where(np.isfinite(flux_all) == True)
#    use_these = np.zeros(length_all)
#    use_these[ifin[0]] = 1

    iuse = np.where(use_these == 1)
    if len(iuse[0]) == 0:
      print('what happened?')

    times = time_all[iuse[0]]
    flux = flux_all[iuse[0]]
    sectors = sector_all[iuse[0]]
    sectors_used = np.unique(sectors).astype(int)     
    sectors_used = list(map(str,sectors_used))


    if (tstar['which_LC'] == 'SAP') & (tstar['SAP_line'] == True):
      slope, intercept = np.polyfit(times, flux, 1)
      sap_line = slope * times + intercept
      flux -= sap_line


#Periodogram Setup
    Pmin = tstar['pmin']
    Pmax = tstar['pmax']
    Fmin = 1/Pmax
    Fmax = 1/Pmin

#	freq_cpm, pow_cpm = LombScargle(lc_cpm["time"], lc_cpm["flux"]).autopower(minimum_frequency=Fmin,maximum_frequency=Fmax)
#	naf= np.array(1/freq_cpm)
#	nap= np.array(pow_cpm)
#	maX = np.argmax(nap)
#	period = (naf[maX])

    time_2 = time.time()
    periods_cpm = np.logspace(np.log10(Pmin),np.log10(Pmax),10000)
    freq_cpm = 1/periods_cpm
    pow_cpm = LombScargle(times, flux).power(freq_cpm)
    period = periods_cpm[np.argmax(pow_cpm)]
    tstar['Prot_LS'] = period
    tstar['Power_LS'] = np.max(pow_cpm)
    
# Amplitude measurement

    perc05 = np.percentile(flux,5)
    perc95 = np.percentile(flux,95)
    amp = float(perc95-perc05)
    tstar['Amplitude'] = amp


# check if double
# read which_LC, then store in that period. 
# store these in star, update Prot_final

    mdub = float(1.0) 
    if tstar['is_it_double'] == 1:
        mdub = float(2.0)
    period_final = float(tstar['Prot_LS']*mdub)
    tstar['Prot_final'] = period_final
    
   
#Figure creation
    if noplot == False:
        panel = plt.figure(constrained_layout=True, figsize= (16,11))
        gs = gridspec.GridSpec(100, 100)

#cpm lightcurve 
# how many sectors?
    #all_sectors = lc_cpm['sector'].unique().astype(int)
#    unique_sectors = sector_all.unique().astype(int)
#    all_sectors = sectors # I'm pretty sure I reran CPM so that this isn't an issue anymore. Bad sectors arent in the CPM file.

    n_sec = len(sectors_used) # this should probably be number used
    n_all_sec = len(unique_sectors)  
    tstar['number_sectors'] = n_sec
    primary_colors = ['b','r']
    color_options = []
    for icol in range(n_sec):
      color_options.append(primary_colors[icol%2])

    n_obs = len(sectors)
    colors = np.repeat('r', n_obs)
    for i in range(n_sec):
        id = np.where(np.array(sectors) == float(sectors_used[i]))
        colors[id] = color_options[i]

    tmin = np.min(times)
    if noplot == False:
        cpmlight = panel.add_subplot(gs[0:40, 0:100])
        cpmlight.set_title(lc_title)
        cpmlight.scatter(times-tmin,flux*100,c = colors,s=15)
        cpmlight.set_xlabel('Day of Observation')
        cpmlight.set_ylabel('Percent Change in Brightness')
        #find midpoint in time array to place amplitude
        amp_time = np.mean(times-tmin)
        #plot amplitude
        cpmlight.plot([amp_time,amp_time],[-amp*100/2,amp*100/2],c='purple')
        if float(tstar['y_max'])>0:
            cpmlight.set_ylim(float(tstar['y_min']),float(tstar['y_max']))
#        if ((float(tstar['x_min']) != 0) | (float(tstar['x_max']) != 0)):
#            cpmlight.set_xlim(float(tstar['x_min'])-0.5,float(tstar['x_max'])+0.5)

        #Mark adding a text for each sector as it is plotted
        bot, top = cpmlight.get_ylim() #find the upper limit for us to put text
        
        for x in sectors.index:
            if x == sectors.index.min():
                cpmlight.text(times[x]-tmin, top*.9, str(int(sectors[x]))) #put sector number for first sector
                cur_sector = sectors[x]
            else:
                if sectors[x]!=cur_sector:
                    cpmlight.text(times[x]-tmin, top*.9, str(int(sectors[x]))) #put sector number for each subsequent sector
                    cur_sector = sectors[x]

# Phased light curve
#cpm_phase = panel.add_subplot(gs[55:, :40])
    if noplot == False:
        cpm_phase = panel.add_subplot(gs[55:, 65:])
        cpm_phase.set_title('Phased Light Curve')
        cpm_phase.scatter(times%period_final,flux*100,c=colors,s=7)
        cpm_phase.set_xlabel('Day in Rotation Cycle')
        #cpm_phase.set_ylabel('Percent Change in Brightness')
        if float(tstar['y_max'])>0:
            cpm_phase.set_ylim(float(tstar['y_min']),float(tstar['y_max']))

#cpm periodogram
    if noplot == False:
        #cpmper = panel.add_subplot(gs[55:,32:60])
        cpmper = panel.add_subplot(gs[55:,34:60])
        cpmper.set_title('Periodogram')
        cpmper.plot(1/freq_cpm, pow_cpm, color = 'black')
        cpmper.set_xlabel('Period (days)')
        cpmper.set_ylabel('Power')
        if tstar['Power_LS']<0.1:
            cpmper.set_yscale('log')
            cpmper.set_ylim(0.001,1)
        if tstar['pxlog'] == 1:
            cpmper.set_xscale('log')
#        cpmper.plot([tstar['Prot_final'],tstar['Prot_final']],[0,1],c='red')
        cpmper.plot(tstar['Prot_final'],0.95,marker='v',markerfacecolor='red',markersize=20,markeredgecolor="black")

#  print(cpmlight.get_xlim())

#    print('First panels: '+str(time.time()-time_3))  	
#   First panels: 0.05 seconds    

#FFI image

    time_4 = time.time()
    if noplot == False:
#        if (tstar['which_sectors'] == 'All'):
        ffi_image = load_ffi_fromfile(tstar['file_ffi'][0])
        if (n_all_sec > 1) & (ffi_test(ffi_image) == 0):
            print('switching sector')
            ffi_image = load_ffi_fromfile(tstar['file_ffi'][1])
        if (tstar['which_sectors'] != 'All') & (np.size(tstar['file_ffi'])>1):
            if tstar['which_sectors'] == '15':
                ffi_image = load_ffi_fromfile(tstar['file_ffi'][0])
            if tstar['which_sectors'] == '16':
                ffi_image = load_ffi_fromfile(tstar['file_ffi'][1])

        ffimage = panel.add_subplot(gs[55:, 0:25])
        ffimage.set_title('TESS Cutout Image')
        color_map = plt.cm.get_cmap('gray')
        reversed_color_map = color_map.reversed()
        ffi_mod = np.clip(ffi_image-np.min(ffi_image),0,1000)
        ffimage.imshow(ffi_mod,origin = 'lower',cmap=reversed_color_map)
        ffimage.plot([15,17],[20,20],color='red')    
        ffimage.plot([23,25],[20,20],color='red')    
        ffimage.plot([20,20],[15,17],color='red')    
        ffimage.plot([20,20],[23,25],color='red')    
        ffimage.set_xlabel('Pixels')
        ffimage.set_ylabel('Pixels')

    if save == True:
#    dir_panels = '/content/gdrive/My Drive/TESS/'+tstar['Cluster']+'/Panels/'
        #dir_panels = '/content/gdrive/My Drive/Plots/Panels_Final/'
        dir_panels = myDir.project_dir(tstar['Cluster'])+'Panels/'
        end = '-User='+tstar['User']+'.png'
        if tstar['Source'][0] == 'G':
            dr2 = tstar['Source'][9:]
        else:
            dr2 = tstar['Source']
        file_save = dir_panels+'GaiaDR2_'+dr2+end
        panel.savefig(file_save, dpi=300, bbox_inches='tight', pad_inches=0.25, transparent=False)

 #   print('Display time:' + str(time.time() - time_1))

    return tstar

####################################################################################
def update_prot_table(table, tstar):
    user = tstar['User']


# (1) where is the star in the table?
    cell = table.find(tstar['Source'])
#	print("Found something at R%sC%s" % (cell.row, cell.col))
    row_number = cell.row

    if user == 'Auto':
        columns = ['Prot','Prot_LS', 'Power_LS', 'TESS_Data']
        n_col = len(columns)
        cols = []
        for column in columns:
            cell = table.find(column)    
            cols.append(cell.col)
        cell_list = [table.cell(row_number,cols[0]),table.cell(row_number,cols[1]),table.cell(row_number,cols[2]),table.cell(row_number,cols[3])]
        cell_list[1].value = tstar['Prot_LS']            
        cell_list[2].value = tstar['Power_LS']            
        if (tstar['exist_ffi'] == 0) or (tstar['exist_cpm'] == 0):
            cell_list[0].value = '-1'
            cell_list[3].value = 'No'
        else:
            cell_list[0].value = ''
            cell_list[3].value = 'Yes'
        table.update_cells(cell_list)

        #added amplitude column, and sector_list 
    if user != 'Auto':
        columns = ['Prot_Final','Prot_LS', 'Power_LS', 'Single_Double', 'Multi', 'Quality', 'LC_Source', 'Class', 'Notes', 'Amp','Sectors_Used','Flares']
        n_col = len(columns)
        cols = [2,3,4,5,6,7,8,9,10,11,12,13]
#        for column in columns:
#            cell = table.find(column)    
#            cols.append(cell.col)
        cell_range = 'B'+str(row_number)+':M'+str(row_number)
        cell_list = table.range(cell_range)
        if tstar['LC_Action'] == 'Publish':
            cell_list[0].value = tstar['Prot_final']
        if tstar['LC_Action'] == 'Good':
            cell_list[0].value = tstar['Prot_final']
        if tstar['LC_Action'] == 'Follow up':
            cell_list[0].value = -tstar['Prot_final']
        if tstar['LC_Action'] == 'Flat':
            cell_list[0].value = 99
        if tstar['LC_Action'] == 'Garbage':
            cell_list[0].value = -99
        cell_list[1].value = tstar['Prot_LS']            
        cell_list[2].value = tstar['Power_LS']            
        cell_list[3].value = tstar['is_it_double']            
        cell_list[4].value = tstar['Multi']            
        cell_list[5].value = tstar['LC_Quality']            
        cell_list[6].value = tstar['which_LC']            
        cell_list[7].value = tstar['LC_Action']                    
        cell_list[8].value = tstar['Notes']
        cell_list[9].value = tstar['Amplitude']
        cell_list[10].value = str(tstar['which_sectors'])
        cell_list[11].value = tstar['Flares']
        table.update_cells(cell_list)

####################################################################################
def stars_todo(table):
    n_stars = len(table)
    not_done = np.zeros(n_stars, dtype=int)
    for i in range(n_stars):
        if len(table['Prot_LS'][i]) == 0:
            not_done[i] = 1
            
    ido = np.where(not_done == 1)
    dr2_list = table['DR2Name'].to_numpy()
    star_list = dr2_list[ido[0]]

    return star_list 

def stars_todo_split(table, user):
    n_stars = len(table)
    not_done = np.zeros(n_stars, dtype=int)
    for i in range(n_stars):
        if len(table['Prot_LS'][i]) == 0:
            not_done[i] = 1
            
    ido = np.where(not_done == 1)
    list = ido[0]
    if user == 'Angeli':
        these = np.where(ido[0] < 202)
        list = ido[0][these]
    if user == 'Isabella':
        these = np.where((ido[0] > 191) & (ido[0] < 379))
        list = ido[0][these]
    if user == 'Linus':
        these = np.where(ido[0] > 373)
        list = ido[0][these]

    dr2_list = table['DR2Name'].to_numpy()
    star_list = dr2_list[list]

    return star_list 

####################################################################################    
def stars_nodata():
    target_table = get_prot_table('Auto')
    target_data = target_table.get_all_records()
    dr2_list = target_table.col_values(1)
    dr2 = np.array(dr2_list[1:])
    num = np.size(dr2)

# load their table
    user_table = get_prot_table('Jason')
# Identify columns        
# Locate Prot_LS column
    cell = user_table.find('Prot_LS')
    col_prot = cell.col
# Locate Power LS column
    cell = user_table.find('Power_LS')
    col_pow = cell.col
# Loop over targets
    for i in range(num):
        row_number = i+2
# Does target have TESS data, according to target table?
        val = target_data[i]['Prot_LS']
        if val == 0:
            user_table.update_cell(row_number,col_prot,0)
            user_table.update_cell(row_number,col_pow,0)
####################################################################################
def stars_nodata_new(Team):
    target_table = get_prot_table('Auto')
    target_data = target_table.get_all_records()
    dr2_list = target_table.col_values(1)
    dr2 = np.array(dr2_list[1:])
    num = np.size(dr2)

# load their table
    user_table = get_prot_table('Jason')
# Identify columns        
# Locate Prot_LS column
    cell = user_table.find('Prot_LS')
    col_prot = cell.col
# Locate Power LS column
    cell = user_table.find('Power_LS')
    col_pow = cell.col
# Loop over targets
    for i in range(num):
        row_number = i+2
# Does target have TESS data, according to target table?
        val = target_data[i]['Prot_LS']
        if val == 0:
            user_table.update_cell(row_number,col_prot,0)
            user_table.update_cell(row_number,col_pow,0)

####################################################################################
def ffi_test(ffi):
    shape = np.shape(ffi)
    val = ffi[int(shape[0]/2),int(shape[1]/2)]
    if np.isfinite(val):
      good_or_bad = 1
    else:
      good_or_bad = 0
    return good_or_bad
####################################################################################
def tess_inspect_not_working(tstar):
    import ipywidgets as widgets
    from ipywidgets import interactive
    from IPython.display import display
    import matplotlib.pylab as pylab
  
    lc_cpm = load_cpm(tstar)
    all_sectors = lc_cpm['sector'].unique().astype(int)
    sectors = sector.unique().astype(int)

    thing_widget = widgets.SelectMultiple(
        options=sectors,
        value=sectors,
        #rows=10,
        description='Sectors',
        disabled=False
        )

    interactive_plot = interactive(idisplay_tess, thing=thing_widget,double=False)
    output = interactive_plot.children[-1]
    interactive_plot
    idisplay()
    return
####################################################################################
####################################################################################
def update_panelname_v1(tstar, locate=False):
    import os
    from glob import glob
    dir_panels = '/content/gdrive/My Drive/Projects/'+tstar['Cluster']+'/Panels/'
    name = dir_panels + '*'+tstar['Source'][9:]+'*'
    file = glob(name)
    if locate == True:
        return np.size(file)
    else:
        end = '-User='+tstar['User']+'-Review='+tstar['LC_Action']+'.png'
        new_file = dir_panels + 'GaiaDR2_'+str(tstar["Source"])[9:]+end
        os.rename(file[0],new_file)
        
def update_panelname(tstar, locate=False):
    import os
    from glob import glob
    dir_panels = myDir.project_dir(tstar['Cluster'])+'Panels/'
    if tstar['Source'][0] == 'G':
        dr2 = tstar['Source'][9:]
    else:
        dr2 = tstar['Source']
    name = dir_panels + '*'+dr2+'*'+tstar['User']+'*'
    file = glob(name)
    if locate == True:
        return np.size(file)
    else:
        end = '-User='+tstar['User']+'-Review='+tstar['LC_Action']+'.png'
        new_file = dir_panels + 'GaiaDR2_'+dr2+end
        os.rename(file[0],new_file)        

#############################################################################
def prot_show_v1(project_name, user, gbr, clusters=False):
#    gbr = target_data['BP_RP'].to_numpy(dtype=float)
    fig1, ax1 = plt.subplots(figsize=(7,6)) 
    ax1.tick_params(axis='both', which='major', labelsize=15)
    aw = 1.5
    ax1.spines['top'].set_linewidth(aw)
    ax1.spines['left'].set_linewidth(aw)
    ax1.spines['right'].set_linewidth(aw)
    ax1.spines['bottom'].set_linewidth(aw)

    prot_table_now = get_prot_table(user,project_name)
    prot_data_val_now = prot_table_now.get_all_values()
    prot_data_now = pd.DataFrame.from_records(prot_data_val_now[1:],columns=prot_data_val_now[0])
    
    pnow = prot_data_now['Prot_Final'].to_numpy()
    uu = np.where((pnow != '') & (pnow != '-1') & (gbr != 'nan'))
    prot_now = pnow[uu[0]]
    color = gbr[uu[0]].astype(float)
    ax1.set_xlim(0.4,2.5)
    ax1.set_xlabel('BP - RP (mag)',fontsize=20)
    ax1.set_ylim(0,20)
    ax1.set_ylabel('Rotation Period (days)',fontsize=20)
    
    if clusters == True:
        file = glob('/content/gdrive/My Drive/Tables/gyro_clusters_draft-2020April08.csv')
        clus = pd.read_csv(file[0])
        indicesPl = np.where((clus["CLUSTER"] ==  "Pleiades") & (clus['BENCH'] == 1))
        indicesPr = np.where((clus["CLUSTER"] == "Praesepe") & (clus['BENCH'] == 1))
        #indicesNGC = np.where((Cluster == "NGC_6811") & (clus['BENCH'] == 1))
        pleiades = clus.iloc[indicesPl]
        praesepe = clus.iloc[indicesPr]
        #NGC6811 = clus.iloc[indicesNGC]
        plt.plot(pleiades["BP_RP"]-0.415*0.12, pleiades["PROT"], markerfacecolor = 'blue', markeredgecolor='black', label = '120 Myr Pleiades',markersize=10,alpha=0.7,linestyle='',marker='.')
        plt.plot(praesepe["BP_RP"]-0.415*0.035, praesepe["PROT"], markerfacecolor = 'cyan', markeredgecolor='black', label = '670 Myr Praesepe',markersize=10,alpha=0.7,linestyle='',marker='.')

    ax1.plot(color, np.array(prot_now,dtype=float),markerfacecolor='red',markeredgecolor='black',marker='*',markersize=20,linestyle='')
    plt.show()
    
    
def prot_show(project_name, user, gbr, clusters=False, pcut=0.0):
#    gbr = target_data['BP_RP'].to_numpy(dtype=float)
    fig1, ax1 = plt.subplots(figsize=(15,9)) 
    ax1.tick_params(axis='both', which='major', labelsize=15)
    aw = 1.5
    ax1.spines['top'].set_linewidth(aw)
    ax1.spines['left'].set_linewidth(aw)
    ax1.spines['right'].set_linewidth(aw)
    ax1.spines['bottom'].set_linewidth(aw)

    prot_table_now = get_prot_table(user,project_name)
    prot_data_val_now = prot_table_now.get_all_values()
    prot_data_now = pd.DataFrame.from_records(prot_data_val_now[1:],columns=prot_data_val_now[0])
    
    pnow = prot_data_now['Prot_Final'].to_numpy()
    qnow = prot_data_now['Quality'].to_numpy()
    cnow = prot_data_now['Class'].to_numpy()

    uu = np.where((pnow != '') & (pnow != '-1') & (gbr != 'nan') & (qnow != '-1') & (cnow == 'Accept'))
#    uu = np.where((pnow != '') & (pnow != '-1') & (gbr != 'nan'))
    prot_now = pnow[uu[0]]
    color = gbr[uu[0]].astype(float)
    power_now = prot_data_now['Power_LS'].to_numpy()
    vv = np.where(power_now[uu[0]].astype(float)>pcut)
    ax1.set_xlim(0.4,3.5)
    ax1.set_xlabel('BP - RP (mag)',fontsize=20)
    ax1.set_ylim(0,25)
    ax1.set_ylabel('Rotation Period (days)',fontsize=20)

    
    if clusters == True:
        file = glob('/content/gdrive/My Drive/Tables/gyro_clusters_draft-2020April08.csv')
        clus = pd.read_csv(file[0])
        indicesPl = np.where((clus["CLUSTER"] ==  "Pleiades") & (clus['BENCH'] == 1))
        indicesPr = np.where((clus["CLUSTER"] == "Praesepe") & (clus['BENCH'] == 1))
        #indicesNGC = np.where((Cluster == "NGC_6811") & (clus['BENCH'] == 1))
        pleiades = clus.iloc[indicesPl]
        praesepe = clus.iloc[indicesPr]
        #NGC6811 = clus.iloc[indicesNGC]
        ax1.plot(pleiades["BP_RP"]-0.415*0.12, pleiades["PROT"], markerfacecolor = 'blue', markeredgecolor='black', label = '120 Myr Pleiades',markersize=10,alpha=0.7,linestyle='',marker='.')
        ax1.plot(praesepe["BP_RP"]-0.415*0.035, praesepe["PROT"], markerfacecolor = 'cyan', markeredgecolor='black', label = '670 Myr Praesepe',markersize=10,alpha=0.7,linestyle='',marker='.')


    ax1.plot(color[vv[0]], np.array(prot_now[vv[0]],dtype=float),markerfacecolor='red',markeredgecolor='black',marker='*',markersize=15,linestyle='')
    print(len(vv[0]))
#    ax1.scatter([1.2215], [11.6770],s=3000,c='green')
#    ax1.plot([1.2375,1.2375],[0,20],c='green')
#    ax1.plot([0.5,2.5],[11.677,11.677],c='green')
    plt.show()
    
#############################################################################
def make_status(project_name, add_user='', add_step = '', change_auto='', reset=False):
    # directory
    dir_project = myDir.project_dir(project_name)
    # ensure the file doesnt already exist
    file_status = glob(dir_project + "Status.txt")

    new_file = 1
    if (np.size(file_status) == 1) | (reset == False):
        bsize = os.path.getsize(file_status[0]) #size in bytes
        if bsize < 40:
            print('remove the file')
            os.remove(file_status[0])  
        else:
            new_file = 0
            status = read_status(project_name)

    if (new_file == 1) | (reset == True):
        status = {'Project':project_name,
            'Users':'Jason Team_Member',
            'Steps':'Status_Initialized',
            'Auto': 'No'}

    if len(add_user) > 0:
        status['Users'] += ' '+add_user

    if len(add_step) > 0:
        status['Steps'] += ' '+add_step

    if len(change_auto) > 0:
        status['Auto'] = change_auto

    lines = []
    # 1: project title
    lines.append('Project: '+project_name+"\n")
    # 2: Users
    lines.append('Users: '+status['Users']+"\n")
    # 3: Steps
    lines.append('Steps: '+status['Steps']+"\n")
    # 4: Has tesscheck_auto been run?
    lines.append('tesscheck_auto: '+status['Auto']+"")
    # 5: Which sectors?
    # 6: Number of repeat sectors?

    # Create the file
    fo = open(dir_project + "Status.txt", "w")
    fo.writelines(lines)
    fo.close()
######################################
def read_status(project_name):
    dir_project = myDir.project_dir(project_name)
    # ensure the file doesnt already exist
    file_status = glob(dir_project + "Status.txt")
    if np.size(file_status) == 0:
        make_status(project_name)
        #print('no file')
        return

    bsize = os.path.getsize(file_status[0]) #size in bytes
    if (bsize < 40):
        os.remove(file_status[0])  
        make_status(project_name)
        #print('no file')
        return

    fo = open(file_status[0], "r")
    lines_current = fo.readlines()
    fo.close()

    # 2 Users
    c_users_line = lines_current[1][0:-1]
    c_users_split = c_users_line.split(sep=' ')
    users_current = c_users_split[1:]

    # 3 Steps
    c_steps_line = lines_current[2][0:-1]
    c_steps_split = c_steps_line.split(sep=' ')
    steps_current = c_steps_split[1:]

    # 4 Auto
    c_line = lines_current[3]
    c_split = c_line.split(sep=' ')
    auto_current = c_split[1:]

    status = {'Project':project_name,
            'Users':list_to_string(users_current),
            'Steps':list_to_string(steps_current),
            'Auto':auto_current[0]}
            
    return status
#############################################################################
def add_user_to_sheet(project_name, user):
    from google.colab import auth
    auth.authenticate_user()
    import gspread
    from oauth2client.client import GoogleCredentials
    gc = gspread.authorize(GoogleCredentials.get_application_default())
    sheet = load_sheet(project_name)
    sheet_list = sheet.worksheets()

    target_table = get_prot_table('Auto',project_name)
    target_data_val = target_table.get_all_values()
    target_data = pd.DataFrame.from_records(target_data_val[1:],columns=target_data_val[0])
    dr2_list = target_data['DR2Name'].to_numpy()
    
    number_of_stars = len(dr2_list)

    sheet_exists = 0
    for isheet in sheet_list:
        if isheet.title == user:
            sheet_exists = 1

    if sheet_exists == 1:
        print('This user has a sheet')
    else:
        print('Making sheet for new user...')
        new_sheet = sheet.add_worksheet(title=user,rows=number_of_stars+1,cols=10)
        columns = ['DR2Name', 'Prot_Final', 'Prot_LS', 'Power_LS', 'Single_Double', 'Multi', 'Quality', 'LC_Source', 'Class', 'Notes','Amp','Sectors_Used','Flares']
        cell_range = 'A1:M1'
        cell_list = new_sheet.range(cell_range)
        i=0
        for cell in cell_list:
            cell.value = columns[i]
            i+=1
        new_sheet.update_cells(cell_list)

        cell_range = 'A2:A'+str(number_of_stars+1)
        cell_list = new_sheet.range(cell_range)
        i=0
        for cell in cell_list:
            cell.value = target_data['DR2Name'][i]
            i+=1
        new_sheet.update_cells(cell_list)

        cell_range = 'B2:B'+str(number_of_stars+1)
        cell_list = new_sheet.range(cell_range)
        i=0
        for cell in cell_list:
            if target_data['Prot'][i] == '-1':
                cell.value = target_data['Prot'][i]
            i+=1
        new_sheet.update_cells(cell_list)

        cell_range = 'C2:C'+str(number_of_stars+1)
        cell_list = new_sheet.range(cell_range)
        i=0
        for cell in cell_list:
            if target_data['Prot_LS'][i] == '0':
                cell.value = target_data['Prot_LS'][i]
            i+=1
        new_sheet.update_cells(cell_list)

        cell_range = 'D2:D'+str(number_of_stars+1)
        cell_list = new_sheet.range(cell_range)
        i=0
        for cell in cell_list:
            if target_data['Power_LS'][i] == '0':
                cell.value = target_data['Power_LS'][i]
            i+=1
        new_sheet.update_cells(cell_list)
    
#    print('Setting up a new user took '+str(time.time()-start_time)+' seconds')


def prot_auto_show(project_name, clusters=False, pcut=0.0, av=0.0):
    fig1, ax1 = plt.subplots(figsize=(15,9)) 
    ax1.tick_params(axis='both', which='major', labelsize=15)
    aw = 1.5
    ax1.spines['top'].set_linewidth(aw)
    ax1.spines['left'].set_linewidth(aw)
    ax1.spines['right'].set_linewidth(aw)
    ax1.spines['bottom'].set_linewidth(aw)

    worksheet = load_sheet(project_name)
    target_table = worksheet.worksheet('Targets')
    target_data_val = target_table.get_all_values()
    target_data = pd.DataFrame.from_records(target_data_val[1:],columns=target_data_val[0])    
    
    gbr = target_data['BP_RP'].to_numpy(dtype=float)
    
    pnow = target_data['Prot_LS'].to_numpy()
    uu = np.where((pnow != '') & (pnow != '-1') & (gbr != 'nan'))
    prot_now = pnow[uu[0]]
    color = gbr[uu[0]].astype(float) - 0.415*av
    power_now = target_data['Power_LS'].to_numpy()
    vv = np.where(power_now[uu[0]].astype(float)>pcut)
    ax1.set_xlim(0.4,3.5)
    ax1.set_xlabel('$(BP - RP)_0$ (mag)',fontsize=20)
    ax1.set_ylim(0,25)
    ax1.set_ylabel('Rotation Period (days)',fontsize=20)

    
    if clusters == True:
        file = glob('/content/gdrive/My Drive/Tables/gyro_clusters_draft-2020April08.csv')
        clus = pd.read_csv(file[0])
        indicesPl = np.where((clus["CLUSTER"] ==  "Pleiades") & (clus['BENCH'] == 1))
        indicesPr = np.where((clus["CLUSTER"] == "Praesepe") & (clus['BENCH'] == 1))
        indicesNGC = np.where((clus["CLUSTER"] == "NGC_6811") & (clus['BENCH'] == 1))
        pleiades = clus.iloc[indicesPl]
        praesepe = clus.iloc[indicesPr]
        NGC6811 = clus.iloc[indicesNGC]
        ax1.plot(pleiades["BP_RP"]-0.415*0.12, pleiades["PROT"], markerfacecolor = 'blue', markeredgecolor='black', label = '120 Myr Pleiades',markersize=10,alpha=0.7,linestyle='',marker='.')
        ax1.plot(praesepe["BP_RP"]-0.415*0.035, praesepe["PROT"], markerfacecolor = 'cyan', markeredgecolor='black', label = '670 Myr Praesepe',markersize=10,alpha=0.7,linestyle='',marker='.')
        ax1.plot(NGC6811["BP_RP"]-0.415*0.15, NGC6811["PROT"], markerfacecolor = 'orange', markeredgecolor='black', label = '1 Gyr NGC 6811',markersize=10,alpha=0.7,linestyle='',marker='.')

    ax1.plot(color[vv[0]], np.array(prot_now[vv[0]],dtype=float),markerfacecolor='lightgreen',markeredgecolor='black',marker='*',markersize=20,linestyle='')
#    ax1.scatter(1.758118-0.415*av,13.2,c='red',s=200)

    print(len(vv[0]))
#    ax1.plot([1.2375,1.2375],[0,20],c='green')
#    ax1.plot([0.5,2.5],[11.677,11.677],c='green')
    plt.show()

