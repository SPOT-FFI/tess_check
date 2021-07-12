"""
Module to contain functions related to the status of a project
"""
import os, sys
from glob import glob

import numpy as np

############# If you move these programs you will need to update these directories and names #############
sys.path.append('/content/gdrive/My Drive/')
from . import myDir as myDir
from .sheet import load_sheet, get_prot_table

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
