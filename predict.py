import sys
import ROOT
import numpy as np
import array
import os
from src.dataloader import *
from catboost import CatBoostClassifier

def predict(rootfile="",
            model_path="",
            tree_name="EventTree"):
    
    # Create a CatBoostClassifier object
    model = CatBoostClassifier()
        
    # Load model from given directory
    model.load_model(model_path)

    # Load MLInput data
    X=load_data(rootfiles=[rootfile], version="predict")

    # No data was successfully loaded
    # We return -1 when data loading fails
    if(type(X)==int):
        print("Failure to load data...")
        return -1

    # Make signal predictions
    prob=model.predict_proba(X)[:,1]
    
    # Load EventTree
    # Here we will create a new weights branch so that the photon-per-photon classification can be stored
    tfile = ROOT.TFile(rootfile, "UPDATE")
    tree = tfile.Get(tree_name)

    # create a new branch for the tree
    weights = array.array('d',100*[0.0])
    weight_branch=tree.Branch("p_gamma", weights,'p_gamma[Nmax]/D')

    # Loop over the events in the EventTree and set the weights array accordingly
    # If the particle is not a photon, set weight to 1
    # If the particle is a photon, set the weight from the prediction
    k=0 # An incrementing variable to step to the next element in "prob" after a photon is found
    N=tree.GetEntries()
    for i in range(N):
        tree.GetEntry(i)
        pid=np.array(tree.pid)
        for j,PID in enumerate(pid):
            if(PID==22):
                weights[j]=prob[k]
                k+=1
            else:
                weights[j]=1
        weight_branch.Fill() # Fill the branch

    # Write the TTree and close the TFile
    tree.Write(tree_name,1) # the "1" forces an overwrite of the previous ttree
    tfile.Close()
        
if __name__ == "__main__":
    args = [arg for arg in sys.argv[1:]]
    predict(*args)
        
        
            
