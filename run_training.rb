#!/usr/bin/env ruby

require 'optparse'

def print_red(text)
  puts "\033[0;31m#{text}\033[0m"
end

def print_green(text)
  puts "\033[0;32m#{text}\033[0m"
end

def create_slurm_file(command, project_directory, job_name)
  slurm_template = <<-SLURM
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=production
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=#{job_name}
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
#SBATCH --output=#{project_directory}/slurm/#{job_name}.out
#SBATCH --error=#{project_directory}/slurm/#{job_name}.err

#{command}
  SLURM

  slurm_file = "#{project_directory}/slurm/#{job_name}.slurm"
  File.open(slurm_file, "w") { |file| file.write(slurm_template) }
  slurm_file
end

def create_gpu_slurm_file(command, project_directory, job_name)
  slurm_template = <<-SLURM
#!/bin/bash
#SBATCH --account=clas12
#SBATCH --partition=scavenger_gpu
#SBATCH --mem-per-cpu=4000
#SBATCH --job-name=#{job_name}
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=#{project_directory}/slurm/#{job_name}.out
#SBATCH --error=#{project_directory}/slurm/#{job_name}.err

#{command}
  SLURM

  slurm_file = "#{project_directory}/slurm/#{job_name}.slurm"
  File.open(slurm_file, "w") { |file| file.write(slurm_template) }
  slurm_file
end

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: run_training.rb [options]"

  opts.on("--project_name PROJECT_NAME", "Specify project name") do |project_name|
    options[:project_name] = "training_projects/"+project_name
  end

  options[:mode] = :normal
  opts.on("--fast", "Run training in fast mode") do
    options[:mode] = :fast
  end

  opts.on("--normal", "Run training in normal mode (default)") do
    options[:mode] = :normal
  end
    
  opts.on("--slurm", "Use SLURM for processing") do
    options[:slurm] = true
  end
end.parse!

if options[:project_name].nil?
  print_red "Error: Project name is missing."
  exit 1
end

project_name = options[:project_name]
project_directory = project_name
save_EventTree = options[:mode] == :fast ? "0" : "1"
use_slurm = options[:slurm]

if use_slurm
  # Create slurm directory if it doesn't exist
  Dir.mkdir("#{project_directory}/slurm") unless Dir.exist?("#{project_directory}/slurm")
end


# Store job IDs of submitted SLURM jobs
job_ids = []

# Read input files from training_hipo_files.txt
count = 0
total = `wc -l < "#{project_directory}/training_hipo_files.txt"`.to_i

File.open("#{project_directory}/training_hipo_files.txt", "r") do |file|
  file.each_line do |hipo_file|
    count += 1

    if use_slurm
      # Create SLURM file for hipo2tree
      slurm_file = create_slurm_file("clas12root -b -q ./hipo2tree.C\\(\\\"#{hipo_file.chomp}\\\",\\\"#{project_directory}/data/#{File.basename(hipo_file.chomp, '.hipo')}.root\\\",#{save_EventTree}\\)", project_directory, "hipo2tree_#{count}")
      # Submit SLURM job for hipo2tree and store job ID
      job_id = `sbatch #{slurm_file}`.strip
      job_ids << job_id.split.last
    else
      # Generate output filename
      outroot = "#{project_directory}/data/#{File.basename(hipo_file.chomp, '.hipo')}.root"
      # Convert hipo to root using clas12root
      system("clas12root -b -q ./hipo2tree.C\\(\\\"#{hipo_file.chomp}\\\",\\\"#{outroot}\\\",#{save_EventTree}\\)")
    end
    
    print_green "Processed #{count}/#{total} files"
  end
end

if use_slurm
  # Create SLURM file for training with dependency on hipo2tree jobs
  slurm_file = create_gpu_slurm_file("/apps/python3/3.9.7/bin/python3 train.py #{project_directory}", project_directory, "train_model")
  dependency = job_ids.empty? ? "" : "--dependency=afterok:#{job_ids.join(',')}"
  # Submit SLURM job for training with dependency on hipo2tree jobs
  system("sbatch #{dependency} #{slurm_file}")
else
  print_green "Training model..."
  system("/apps/python3/3.9.7/bin/python3 train.py #{project_directory}")
end

if options[:mode] == :normal
  print_green "Making predictions..."
  Dir.glob("#{project_directory}/data/*").each do |file|
    print_green "Predicting for file: #{file}"
    system("/apps/python3/3.9.7/bin/python3 predict.py #{file} #{project_directory}/model/catboost_model")
  end

  print_green "Processing Diphotons..."

  Dir.glob("#{project_directory}/data/*").each do |file|
    print_green "Processing Diphotons for file: #{file}"
    system("clas12root -b -q ./buildDiphotons.C\\(\\\"#{file}\\\"\\)")
  end
end

print_green "Finished processing all files."
