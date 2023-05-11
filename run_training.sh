#!/bin/bash

print_red () {
  echo -e "\033[0;31m$1\033[0m"
}

print_green () {
  echo -e "\033[0;32m$1\033[0m"
}

if [ -z "$1" ]; then
    print_red "Error: project name is missing"
    print_red "Usage: ./run_training.sh [project_name]"
    print_red "\n\tProjects are created by the ./create_training_project.sh script"
    exit 1
fi

project_name="$1"
project_directory="training_projects/$project_name"

# Read input files from training_hipo_files.txt
count=0
total=$(wc -l < "$project_directory/training_hipo_files.txt")
while read hipo_file; do
    count=$((count+1))
    # Generate output filename
    outroot="$project_directory/data/$(basename "$hipo_file" .hipo).root"
    
    # Convert hipo to root using clas12root
    clas12root -b -q ./hipo2tree.C\(\"${hipo_file}\",\"${outroot}\"\)
    clas12root -b -q ./EventTree2MLinput.C\(\"${outroot}\"\)
    
    print_green "Processed $count/$total files"
done < "$project_directory/training_hipo_files.txt"

print_green "Training model..."
/apps/python3/3.9.7/bin/python3 train.py $project_directory

print_green "Making predictions..."
for file in "$project_directory/data"/*; do
    print_green "Predicting for file: $file"
    /apps/python3/3.9.7/bin/python3 predict.py "$file" "$project_directory/model/catboost_model"
done

print_green "Processing Diphotons..."

for file in "$project_directory/data"/*; do
    print_green "Processing Diphotons for file: $file"
    clas12root -b -q ./buildDiphotons.C\(\"$file\"\)
done

print_green "Finished processing all files."
