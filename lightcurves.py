"""
Module to contain functions related to SAP and CPM light curves
"""
import sys
from glob import glob

import pandas as pd

############# If you move these programs you will need to update these directories and names #############
sys.path.append('/content/gdrive/My Drive/')
from .import myDir as myDir


####################################################################################
# Functions to locate a star's TESS data

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
