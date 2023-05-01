# CLAS12 Photon Classification
This repository contains scripts and code to perform photon classification on data files from the CLAS12 experiment at Jefferson Lab.

Pretrained Gradient Boosted Trees (see https://catboost.ai/en/docs/) models are provided using 10.6 GeV inbending (45nA bkg), and 10.6 GeV outbending (50nA bkg) clasdis Monte Carlo. Users can also train their own models, which may need to be done if the default cuts included in this repository are loosened/tweaked. 

## Tutorial
---

### Requirements

Ensure you have `clas12root` installed on ifarm and have the capability to run `module load clas12/pro`

The following python packages must be installed under python3.9.7. To prepare this, do the following...

  1. module unload python
  2. module load python3/3.9.7
  3. pip3 install copy, pyroot, catboost, os, numpy, array, sys, uproot4, pandas, scikit-learn, matplotlib, seaborn, tqdm
  
### Predicting with pre-trained models

The quickest implementation of this repository is to run one of the pre-trained RG-A photon classifiers on a RG-A experimental data `.hipo` file. The repository includes a script called `classify_hipo.sh`, which automates this pipeline for a single given .hipo file and selected model. Symbolic links to the RG-A experimental data are provided in the repository. As an example:

```
./classify_hipo.sh
```

will print out the general usage. If the requirements are met, the following one-liner will analyze a `.hipo` file with the photon classifier:

```
./classify_hipo.sh fall2018_rga_inbending_nSidis/nSidis_005032.hipo pretrained_models/model_rga_inbending outroot
```

The output of `./classify_hipo.sh` will be a `.root` file in `outroot/` containing three TTrees. These are *EventTree*, *MLinput*, and *diphoton*. See below for more info on what is stored in each TTree. The *diphoton* TTree contains event information on a diphoton-by-diphoton level, such as the transverse momentum of the diphoton, its energy, the event x,Q2,W, etc. and importantly, the classifier output for each of the photons. One could plot the `M_gg` distribution for the `diphoton` TTree with a classifier threshold of p>0.9 as follows:

  ```
  diphoton->Draw("M_gg" , "p_gamma_1 > 0.9 && p_gamma_2 > 0.9") 
  ```

### Training New Models

The first step to take when training a new model is to run `./create_training_model.sh`. This script will prompt the user for a project name, and will create a subdirectory in `./training_models/<PROJECT_NAME>/`. Two further subdirectories are created, one storing the "data" (i.e. root files) and the other storing model-specific information (roc curves, confusion matrices, feature importances). 

The second step is to edit the `model_params.yaml` and `training_hipo_files.txt` in the new project directory, following the instructions outputted by `./create_training_model.sh`. 

The third step, which could be time intensive depending on the number of Monte Carlo files used for training, is to run `./run_training.sh [PROJECT_NAME]`. This program will pass the Monte Carlo hipo files through the pipeline (discussed in detail below), adding a training step through `train.py`. The trained model will be saved in `./trained_models/<PROJECT_NAME>/model/catboost_model` and can be used as the second argument in `./classify_hipo.sh` when wanting to test it on experimental data. 

## Source Code and Macros
---

The pipeline involves several steps:

1. The script `hipo2tree.C` is run on a `.hipo file` to create a `.root` file in `outroot/` containing a TTree called `EventTree`. This TTree stores branches that are either single valued (x,Q2,W,etc.) or arrays of size Nmax, where Nmax is the number of particles in an event and thus the length of the arrays (px, theta, E, etc.). 
2. Next, `EventTree2MLinput.C` is run on the `.root` file to create a new TTree called `MLinput`. The `MLinput` TTree contains the input parameters for the pre-trained Gradient Boosted Trees model, which will classify each photon.`
3. The `predict.py` script is then run on the `MLinput` TTree using a selected model to perform photon classification, and the results are stored in a new array branch called `p_gamma` in the original `EventTree`. The values for `p_gamma` range from 0 to 1 for photons, where 1 is the maximum likelihood (according to the model prediction) that the photon is signal. Since `p_gamma` is an array, where each element represents a different particle in the event, `p_gamma` elements are naturally set to 1 for particles that are not identified by the EventBuilder as `pid==22`.
4. Finally, the `buildDiphotons.C` script is run on the `EventTree` to create a new TTree called `diphoton`, whose branches are only single-valued (x,Q2,z,pT,M_gg, etc). The key branches here are `p_gamma_1` and `p_gamma_2`, which represents the classifier output for each photon of the diphoton, respectively.

There are currently two trained GBT models provided in the repository, located under `trained_models/`. These are "model_rga_inbending" and "model_rga_outbending". They are pre-trained using Monte Carlo simulations at `/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus-1/v1/bkg45nA_10604MeV/` and `/cache/clas12/rg-a/production/montecarlo/clasdis/fall2018/torus+1/v1/bkg50nA_10604MeV/`. 

Additional code is stored in `src/` to facilitate applying cuts (`CutManager`) , establishing a particle structure containing REC::Particle and REC::Calorimeter info (`Structs.h`) , and parsing through the hipo banks to fill these structs (`HipoBankInterface`).

Lastly, an example jupyter notebook is provided to showcase the performance of one of the pretrained classifiers.




