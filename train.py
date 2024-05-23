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
    
    # Path to the save CatBoostModel
    model_path = projectdir+"/model/catboost_model"
    
    # Parameters for the model (ex: tree depth)
    model_params = load_params(yamlfile)
    
    
    # Create CatBoost Model
    model = CatBoostClassifier(**model_params,
                    custom_loss=[metrics.Accuracy()], 
                    random_seed=42,
                    #task_type="GPU",
                    #devices='0:1')
                    task_type="CPU")
    
    # Load the rootfiles for the learning
    rootfiles = load_files(projectdir)
    print(len(rootfiles),"found for training")
    
    # Load in data into a training and validation sample
    Nfiles = len(rootfiles)
    for iroot,rootfile in enumerate(rootfiles):
        

        
        print("Loading data for file <"+rootfile+"> ("+str(iroot+1)+"/"+str(Nfiles)+")")
        try:
            train_pool, validation_pool = load_data(rootfiles = [rootfile],
                                                    version="train",
                                                    split_ratio = 0.75,
                                                    random_seed = 42)
        except:
            continue
        print("\t Starting training for file <"+rootfile+">")
        X_train = train_pool[0]
        y_train = train_pool[1]
        X_validation = validation_pool[0]
        y_validation = validation_pool[1]
        numeric_train_pool = Pool(X_train, y_train)
        numeric_val_pool = Pool(X_validation, y_validation)
        
        # Perform the training
        if iroot==0:
            model.fit(numeric_train_pool, verbose=1,eval_set=numeric_val_pool,early_stopping_rounds=50)
        else:
            model.fit(numeric_train_pool, verbose=1,eval_set=numeric_val_pool,early_stopping_rounds=50,init_model=model_path)
                
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