import pandas as pd
import numpy as np
import uproot4 as uproot
import os

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

# Inputs:
#   rootfiles = [file1.root, file2.root, ...] created by photonML.C
#   ttree     = "MLInput" or any other TTree name if needed
#   version   = Determines if data must be split for training
#   split_ratio = Ratio for splitting data for training and testing
#   random_seed = Seed for random number generator
def load_data(rootfiles=[""]):

    if(type(rootfiles)!=list):
        rootfiles=[rootfiles]
        print("WARNING: Need to convert <rootfiles> to list")
    
    df=pd.DataFrame()

    ttree="MLinput"
    
    for ifile,rootfile in enumerate(rootfiles):

        # Open the file in uproot
        u = uproot.open(rootfile)
        found=False
        for key in u.keys():
            if(ttree in key):
                found=True
                break
        
        if(not found):
            print("Skipping file",rootfile,"...missing TTree=",ttree)
            continue
            
        try:
            u = u["MLinput"]
        except:
            print("Skipping file",rootfile,"...missing TTree=",ttree)
            continue
        # If the TTree is empty, continue
        if(u.num_entries==0):
            continue

        # Get dataframe params
        
        branchnames,keys,nns = get_keys(u)
        if(ifile==0): 
            # Create dataframe from first file
            # This dataframe will be appended to for each file
            df = pd.DataFrame(columns=keys)
        
        # Create temporary dataframe for each file
        tmp_df = pd.DataFrame(columns=keys)
        
        # Fill temporary dataframe
        for b,k,n in zip(branchnames,keys,nns):
            if(n==-1):
                tmp_df[k]=u[b].array(library="np")
            else:
                tmp_df[k]=np.array(u[b].array(library="ak")[:,n],dtype=float)
                
        # Concatenate with main dataframe
        df=pd.concat([df,tmp_df], ignore_index=True,axis=0)

        # Delete temporary dataframe
        del tmp_df

    # If the dataframe is empty, return -1
    if(df.empty):
        return -1

    # Create dataset for evaluation
    return df
