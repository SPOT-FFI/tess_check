"""
Module to contain functions related to opening and accessing spreadsheets
"""
import os, sys

import numpy as np
import gspread
from google.colab import auth
from oauth2client.client import GoogleCredentials

############# If you move these programs you will need to update these directories and names #############
sys.path.append('/content/gdrive/My Drive/')
from . import myDir as myDir



####################################################################################
# Load sheet
def load_sheet(project_name):
    auth.authenticate_user()
    gc = gspread.authorize(GoogleCredentials.get_application_default())

    dir_project = myDir.project_dir(project_name)
    sheet_id_file = os.path.join(dir_project,f"{project_name}.txt")

    f = open(sheet_id_file,"r")
    doc_id = f.readline()
    f.close()
    sheet = gc.open_by_key(doc_id)
    return sheet


####################################################################################
def get_prot_table(user,project_name):
	sheet_name = project_name
	auth.authenticate_user()
	gc = gspread.authorize(GoogleCredentials.get_application_default())
	# load the table
	worksheet = load_sheet(sheet_name)
	if user == 'Auto':
		table = worksheet.worksheet('Targets')
	else:
		table = worksheet.worksheet(user)
	return table

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

    dr2_list = table['DR2Name'].to_numpy()
    star_list = dr2_list[list]

    return star_list

####################################################################################
# TODO: Why is Jason a hard-coded user here? Are we still using this function?
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
# TODO: Why is Jason a hard-coded user here? Are we still using this function?
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
