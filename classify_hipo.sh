#!/bin/bash
# Define horizontal line
hl="================================================================================"

# Define ANSI escape codes for colors
GREEN='\033[0;32m'
RED='\033[0;31m'
NC='\033[0m' # No Color

# Define print_green function
print_green() {
  echo -e "${GREEN}$1${NC}"
}

# Define print_red function
print_red() {
  echo -e "${RED}$1${NC}"
}



# Check if the current directory's basename is clas12_photon_classifier
if [ "$(basename $(pwd))" != "clas12_photon_classifier" ]; then
  print_red "Error: This script must be run from the clas12_photon_classifier directory"
  exit 1
fi

# Check if the number of arguments is less than 2
if [ "$#" -lt 2 ]; then
  # Print usage statement and available model options
  echo "Usage: "
  echo "  ./classify_hipo.sh [.hipo file] [model_name]"
  echo ""
  echo "Description: "
  echo "  This script takes a hipo file as input and passes it through a CLAS12ROOT+PYROOT pipeline to perform photon classification."
  echo "  [.hipo file] --> hipo2tree.C --> [.root file w/ EventTree] --> EventTree2MLinput.C --> [.root file w/ MLinput TTree] --> predict.py --> [p_gamma classifier branch in EventTree] --> buildDiphotons.C --> [diphoton TTree created in .root file]"
  echo ""
  echo "Arguments: "
  echo "  [.hipo file] --> absolute/relative path to a single .hipo file (ex: fall2018_rga_inbending_nSidis/nSidis_005032.hipo)"
  echo "  [model_name] --> Which model to use for photon classification. Available options are:"
  for model_path in trained_models/*
  do
    echo "                   - $(basename $model_path)"
  done
  # Exit the script
  exit
fi


# Get the input file and model paths
hipo_file=$1
model=trained_models/$2

# Check if the input .hipo file exists
if ! test -f "$hipo_file"; then
  print_red "Error: Input .hipo file does not exist"
  exit 1
fi

# Check if the selected model exists
if ! test -f "$model"; then
  print_red "Error: Selected model does not exist"
  exit 1
fi

# Print success message in green if input .hipo file and selected model exist
print_green "Input .hipo file and selected model exist"

# Create the output .root file path
outroot=outroot/$(basename "$hipo_file" .hipo).root

################################################################################
# Run hipo2tree.C on hipo file
################################################################################

echo $hl
print_green "Running hipo2tree.C on $(basename "$hipo_file")"
clas12root -b -q ./hipo2tree.C\(\"${hipo_file}\",\"${outroot}\"\)
print_green "Root file $outroot created containing TTree \"EventTree\")"


################################################################################
# Run EventTree2MLinput.C on root file
################################################################################

echo $hl
print_green "Running EventTree2MLinput.C on $outroot)"
clas12root -b -q ./EventTree2MLinput.C\(\"${outroot}\"\)
print_green "New TTree \"MLinput\" created in $outroot"

################################################################################
# Read in the MLinput TTree to perform photon classification
################################################################################

echo $hl
print_green "Running predict.py on $outroot using model $model"
echo -e "\t Unloading modules root, python"
source /etc/profile.d/modules.sh
module unload root
module unload python
echo -e "\t Loading python3/3.9.7 and latest root module"
module load python3/3.9.7
module load root
/apps/python3/3.9.7/bin/python3 predict.py $outroot $model
print_green "Prediction done. New photon classification branch \"p_gamma\" stored in EventTree for $outroot"
echo -e "Unloading root module"
module unload root
echo -e "Reloading clas12/pro"
module load clas12/pro


################################################################################
# Create the diphotons from the EventTree
################################################################################
echo $hl
print_green "Running buildDiphotons.C on $outroot"
clas12root -b -q ./buildDiphotons.C\(\"${outroot}\"\)
print_green "New TTree \"diphotons\" created in $outroot"
echo $hl

echo -e "\n\n All done \n\n"