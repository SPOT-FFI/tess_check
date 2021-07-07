#import tess_cpm
import numpy as np
from glob import glob
from astroquery.vizier import Vizier
from astroquery.mast import Tesscut
from astropy.coordinates import SkyCoord
import astropy.units as u
#import matplotlib.pyplot as plt


# this is the Gaia DR2 cluster membership catalog
def cluster_catalog(cluster_name):
	member_catalog_code = "J/A+A/618/A93/members" 
	quereycatalog = Vizier(column_filters = {"Cluster":cluster_name})
	quereycatalog.ROW_LIMIT = -1
	cluster = quereycatalog.get_catalogs(member_catalog_code)
	return cluster


# this is the Gaia DR2 cluster membership catalog
#def cluster_download(cluster_name):
##	tdir = '/Users/srmp/Dropbox/TESS/'+cluster_name+'/'
#	tdir = '/Users/jasoncurtis/Dropbox/TESS/'+cluster_name+'/'

#	member_catalog_code = "J/A+A/618/A93/members" 
#	quereycatalog = Vizier(column_filters = {"Cluster":cluster_name})
#	quereycatalog.ROW_LIMIT = -1
#	cluster = quereycatalog.get_catalogs(member_catalog_code)
#	number_of_stars = len(cluster[0])
#	RA = cluster[0]['RA_ICRS']
#	Dec = cluster[0]['DE_ICRS']
#	for i in range (0,number_of_stars):
#		cutout_coord = SkyCoord(RA[i], Dec[i],unit = "deg")
#		manifest = Tesscut.download_cutouts(cutout_coord, size=40, path=tdir)
#		print(i,number_of_stars)


# Run CPM tools