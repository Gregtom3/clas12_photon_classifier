import sys
import os
from src.dataloader import *
from src.handler import *
from src.plotter import *
from catboost import CatBoostClassifier, Pool, metrics, cv

def train(projectdir = ""):
    
    # Load the parameters from the yamlfile
    yamlfile=projectdir+"/model_params.yaml"
    if(not os.path.exists(yamlfile)):
        print("YAML file",yamlfile,"not found...Aborting...")
        return -1
   
    model_params = load_params(yamlfile)
    
    # Load the rootfiles for the learning
    rootfiles = load_files(projectdir)
    print(len(rootfiles),"found for training")
    
    # Load in data into a training and validation sample
    print("Loading data for",len(rootfiles),"root files")
    train_pool, validation_pool = load_data(rootfiles = rootfiles,
                                            version="train",
                                            split_ratio = 0.75,
                                            random_seed = 42)
    

    X_train = train_pool[0]
    y_train = train_pool[1]
    X_validation = validation_pool[0]
    y_validation = validation_pool[1]
    
    numeric_train_pool = Pool(X_train, y_train)
    numeric_val_pool = Pool(X_validation, y_validation)
    
    # Create CatBoost Model
    model = CatBoostClassifier(**model_params,
                            custom_loss=[metrics.Accuracy()], 
                            random_seed=42,
                            #task_type="CPU"
                            task_type="CPU")
                            #devices='0:1')

    # Perform the training
    print("Starting training")
    model.fit(numeric_train_pool, verbose=1,eval_set=numeric_val_pool)
    
    model.save_model(projectdir+"/model/catboost_model")

    #Import the model and perform predictions with the validation set
    trained_model = CatBoostClassifier()
    trained_model.load_model(projectdir+"/model/catboost_model")

    # Make predictions from the model
    X_validation = validation_pool[0]
    y_validation = validation_pool[1]
    predictions = make_predictions(trained_model,X_validation)

    # Create plots and save them to the model directories
    make_plots(X_validation,y_validation,predictions,projectdir)

    # Save the parameter importances
    feature_names=X_validation.keys()
    save_feature_importance(trained_model,feature_names,projectdir)

if __name__ == "__main__":
    args = [arg for arg in sys.argv[1:]]
    train(*args)