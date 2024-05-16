import pandas as pd
import numpy as np
import uproot4 as uproot
import yaml
import os
from sklearn.model_selection import train_test_split
from tqdm import tqdm

# This function processes a dictionary to determine the keys, original keys, and neighbor parameters of the dictionary. It returns a list of keys, a list of original keys, and a list of neighbor parameters.
# Input: uproot TTree
# Outputs: List of key names, list of modified key names, and index of neighbor
def get_keys(u):
    
    m_g = u["m_g"].array(library="np")[0]
    m_ch = u["m_ch"].array(library="np")[0]
    m_nh = u["m_nh"].array(library="np")[0]
    
    keys = u.keys()
    keys = [key for key in keys if not key in ["m_g","m_ch","m_nh"]]
    
    original_keys=[]
    modified_keys=[]
    nns = []
    
    for key in keys:
        M=-1
        # Since the model can train on multiple neighbor parameters, ensure there
        # is a unique input for each neighbor
        if("_gamma" in key):
            M=m_g
        elif("_ch" in key):
            M=m_ch
        elif("_nh" in key):
            M=m_nh
            
        if(M==-1):
            original_keys.append(key)
            modified_keys.append(key)
            nns.append(-1)
        else:
            for i in range(M):
                original_keys.append(key)
                modified_keys.append(key+"_"+str(i))
                nns.append(i)
    
    return original_keys,modified_keys,nns

# Creates dataset for training or evaluation
# This function reads in root files, extracts the TTree and converts it into a dataframe, then splits the data into a training and validation set according to the given split ratio and random seed. The output is four objects, X_train, X_validation, y_train, y_validation containing the training and validation features and labels, respectively.

def load_data(rootfiles=[""], version="train", split_ratio=0.75, random_seed=42):
    ttree = "MLinput" 
    assert(version in ["train", "predict"])
    
    if not isinstance(rootfiles, list):
        rootfiles = [rootfiles]
        print("WARNING: Need to convert <rootfiles> to list")
    
    data_list = []
    
    # Loop over rootfiles
    for ifile, rootfile in enumerate(rootfiles):
        # Open the file in uproot
        try:
            u = uproot.open(rootfile)
            tree = u[ttree]
        except KeyError:
            print(f"Skipping file {rootfile}...missing TTree={ttree}")
            continue
        except:
            print(f"Skipping file {rootfile}...unexpected error with TTree={ttree}")
            continue

        if tree.num_entries == 0:
            continue
        
        # Get dataframe params
        branchnames, keys, nns = get_keys(tree)
        arrays = tree.arrays(branchnames, library="ak")
        # Convert awkward arrays to pandas DataFrame directly
        tmp_df = pd.DataFrame({k: (np.array(arrays[b]) if n == -1 else np.array(arrays[b][:, n], dtype=np.float32)) for b, k, n in zip(branchnames, keys, nns)})
        data_list.append(tmp_df)
    # Concatenate all DataFrames in the list
    if not data_list:
        return -1
    
    df = pd.concat(data_list, ignore_index=True)
    # Create dataset for training/evaluation
    X = df.drop("photon_has_match", axis=1)
    
    if version == "predict":
        return X
    else:
        y = df["photon_has_match"]
        X_train, X_validation, y_train, y_validation = train_test_split(X, y, train_size=split_ratio, random_state=random_seed)
        return [X_train, y_train], [X_validation, y_validation]
    
    
    
#This code opens a yaml file, and creates a dictionary of the parameters stored in the file.
def load_params(inyaml=""):
    with open(inyaml, 'r') as yaml_file:
        params = yaml.load(yaml_file, Loader=yaml.FullLoader)

    return params

#This function is used to load the files for the machine learning algorithm.
def load_files(projectdir=""):

    rootdir=projectdir+"/data/"
    
    root_files = []

    for file in os.listdir(rootdir):
        if (file.endswith(".root")):
            root_files.append(rootdir+"/"+file)
    
    return root_files
