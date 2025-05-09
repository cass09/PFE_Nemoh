# Imports
from math import pi
import numpy as np
from glob import glob
import os
import shutil


def nemoh_structure(folder):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    else:
        results = os.path.join(folder, "results")
        mesh_ = os.path.join(folder, "mesh")

        if folder_empty(results)==False:
            override = input("Results directory is not empty: override [y]/n?")
            if override.lower() != 'n':
                folder_clear(results)
        if folder_empty(results)==True:
            os.rmdir(results)
        
        if folder_empty(mesh_)==False:
            override = input("Mesh directory is not empty: override [y]/n?")
            if override.lower() != 'n':
                folder_clear(mesh_)
        if folder_empty(mesh_)==True:
            os.rmdir(mesh_)

def file_exist(fn):
    return os.path.isfile(fn)

def folder_empty(folder):
    if os.path.isdir(folder) and os.path.exists(folder):
        if len(os.listdir(folder)) == 0:
            return True
        else:    
            return False                
    else:
        return -1

def folder_clear(folder):
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))
