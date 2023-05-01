#!/bin/bash
echo -e "\n\n\n"
echo "Enter the name of your project: "
read project_name

project_directory="training_projects/$project_name"

if [ -d "$project_directory" ]; then
    echo "Directory already exists. Do you want to overwrite it? (y/n)"
    read confirm_overwrite
    if [ "$confirm_overwrite" = "y" ]; then
        rm -rf "$project_directory"
        mkdir "$project_directory"
        mkdir "$project_directory/data"
        mkdir "$project_directory/model"
        echo "Directory overwritten"
    else
        echo "Operation cancelled"
    fi
else
    mkdir "$project_directory"
    mkdir "$project_directory/data"
    mkdir "$project_directory/model"
    echo "Directory created"
fi

cat << EOF > $project_directory/model_params.yaml
learning_rate: 0.1
depth: 5
iterations: 100
min_data_in_leaf: 5
EOF

cat << EOF > $project_directory/training_hipo_files.txt
mc_rga_inbending/45nA_job_3051_0.hipo
mc_rga_inbending/45nA_job_3051_1.hipo
EOF
echo -e "========================================================================\n"
echo "TO DO"
echo -e "\n========================================================================\n"
echo -e "1. Open $project_directory and edit the \"training_hipo_files.txt\".  Write the list of Monte Carlo files for training. The path to the Monte Carlo files can be absolute OR relative starting from the main repository directory (using the included symlinks)\n"

echo -e "2. Within $project_directory edit the \"model_parms.yaml\" file to manually edit the CatBoost training model parameters (see https://catboost.ai/en/docs/references/training-parameters/)\n"

echo -e "3. After the above two are edited, run [./run_training.sh $project_name]"
echo -e "\n========================================================================\n"
