"""
Module to contain functions related to spreadsheet users
"""

import os, sys

import numpy as np
import pandas as pd
import gspread
from google.colab import auth
from oauth2client.client import GoogleCredentials

############# If you move these programs you will need to update these directories and names #############
sys.path.append('/content/gdrive/My Drive/')
from . import myDir as myDir
from .status import read_status, make_status
from .sheet import load_sheet, get_prot_table

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


#############################################################################
def add_user_to_sheet(project_name, user):
    auth.authenticate_user()
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
