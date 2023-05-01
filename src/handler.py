import os
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier

#This function makes predictions
def make_predictions(model=0, x=[]):
    return np.array(model.predict_proba(x)[:,1],dtype=float)
    
# This function saves the feature importances of a given model
def save_feature_importance(model=0, feature_names=[], projectdir=""):
    # Get feature importances
    importances = model.get_feature_importance()
    # Save feature importance to dataframe
    dfPars = pd.DataFrame(data={"Parameter": feature_names,"Importance": importances})
    # Sort by most important
    dfPars = dfPars.sort_values(by="Importance",ascending=False)
    # Save to csv
    dfPars.to_csv(projectdir+"/model/param_importance.csv",index=False)
    
