import os, shutil

def nemoh_structure(folder):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    else:
        results = os.path.join(folder, "results")
        Motion = os.path.join(folder, "Motion")
        resultsIT = os.path.join(folder, "resultsIT")
        mesh_ = os.path.join(folder, "mesh")
        Operators = os.path.join(folder, "Operators")
        if folder_empty(results)==False:
            override = input("Results directory is not empty: override [y]/n?")
            if override.lower() != 'n':
                folder_clear(results)
                folder_clear(resultsIT)
                folder_clear(Motion)
                folder_clear(Operators)
        if folder_empty(results)==True:
            os.rmdir(results)
            os.rmdir(resultsIT)
            os.rmdir(Operators)
            os.rmdir(Motion)
        if folder_empty(mesh_)==False:
            override = input("Mesh directory is not empty: override [y]/n?")
            if override.lower() != 'n':
                folder_clear(mesh_)
        if folder_empty(mesh_)==True:
            os.rmdir(mesh_)
        for f in os.listdir(folder):
            if f.endswith(".p"):
                try:
                    os.remove(os.path.join(folder, f))
                except OSError as e:
                    print(f"Problem {f} : {e}")
def nemoh_clean(folder):
    
    results = os.path.join(folder, "results")
    mesh_ = os.path.join(folder, "mesh")
    # shutil.rmtree(os.path.join(folder, results))
    if os.path.isdir(mesh_):
        shutil.rmtree(os.path.join(folder, mesh_))
    for f in os.listdir(folder):
        if f.endswith("Normalvelocities.dat"):
            try:
                os.remove(os.path.join(folder, f))
            except OSError as e:
                print(f"Problem {f} : {e}")

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
