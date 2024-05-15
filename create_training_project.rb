#!/usr/bin/env ruby
require 'fileutils'
require 'time'

puts "Select a run group (rga, rgb, rgc):"
run_group = gets.chomp.downcase

# List directories and filter by run group
dirs = Dir.entries('.').select { |entry| File.directory?(entry) && entry.include?(run_group) }

puts "Select bending option (inbending, outbending):"
bending = gets.chomp.downcase

puts "Select pass (pass1, pass2):"
pass = gets.chomp.downcase



# Create an array with the current configuration
current_config = [pass, run_group, bending]

'pass1_mc_rga_inbending'
'pass1_mc_rga_outbending'
'pass2_mc_rga'
'pass1_mc_rgc'

# Define the path based on combinations
path = case current_config
when ['pass1', 'rga', 'inbending']
    'pass1_mc_rga_inbending'
when ['pass1', 'rga', 'outbending']
    'pass1_mc_rga_outbending'
when ['pass2', 'rga', 'inbending']
    'pass2_mc_rga'
when ['pass2', 'rga', 'outbending']
    'pass2_mc_rga'
when ['pass1', 'rgb', 'inbending']
    'pass1_mc_rga_inbending'
when ['pass1', 'rgb', 'outbending']
    'pass1_mc_rga_outbending'
when ['pass2', 'rgb', 'inbending']
    'pass2_mc_rga'
when ['pass2', 'rgb', 'outbending']
    'pass2_mc_rga'
when ['pass1', 'rgc', 'inbending']
    'pass1_mc_rgc'
when ['pass1', 'rgc', 'outbending']
    'pass1_mc_rgc'
when ['pass2', 'rgc', 'inbending']
    'pass2_mc_rga'
when ['pass2', 'rgc', 'outbending']
    'pass2_mc_rga'
else
    # Default case if no matches
    puts "No matching configuration for #{current_config.inspect}."
    return
end


# Adjust file pattern for special cases
file_pattern = "*.hipo"
if path == 'pass2_mc_rga'
  file_pattern = bending == 'inbending' ? "inb-*.hipo" : "outb-*.hipo"
end

# Count .hipo files according to the adjusted pattern
hipo_files = Dir.glob("#{path}/#{file_pattern}")
puts "Found #{hipo_files.size} .hipo files in #{path}"

puts "How many files do you want to analyze?"
num_files = gets.chomp.to_i

# Prepare the base of the project directory
base_subdirectory = "#{run_group}_#{bending}_#{pass}_version"
base_directory = "training_projects/#{base_subdirectory}"
version = 0

# Construct the initial project directory name
timestamp = Time.now.strftime("%Y%m%d%H%M%S")
project_directory = "#{base_directory}#{version}_#{timestamp}"

# Find the highest version number in directories starting with the base_directory
existing_versions = Dir.glob("#{base_directory}*").map do |dir|
  dir.match(/version(\d+)_/)[1].to_i if dir.match(/version(\d+)_/)
end.compact

# If directories already exist, pick the next version number
unless existing_versions.empty?
  version = existing_versions.max + 1
  project_directory = "#{base_directory}#{version}_#{timestamp}"
end

# Create the project directory
FileUtils.mkdir_p(project_directory)
FileUtils.mkdir_p("#{project_directory}/data")
FileUtils.mkdir_p("#{project_directory}/model")
# Write the selected hipo files to a configuration file
File.open("#{project_directory}/training_hipo_files.txt", 'w') do |file|
    hipo_files.first(num_files).each { |file_path| file.puts file_path }
end

File.open("#{project_directory}/model_params.yaml", "w") do |file|
  file.puts "learning_rate: 0.1"
  file.puts "depth: 5"
  file.puts "iterations: 100"
  file.puts "min_data_in_leaf: 5"
end

puts "Project setup complete. Directory: #{project_directory}"


puts "========================================================================\n"
puts "TO DO\n"
puts "\n========================================================================\n"
puts "1. Within #{project_directory} edit the \"model_params.yaml\" file to manually edit the CatBoost training model parameters (see https://catboost.ai/en/docs/references/training-parameters/)\n"
puts "2. After the above two are edited, run [./run_training.rb --project_name #{base_subdirectory}#{version}_#{timestamp}]"
puts "    Preview options with --help."
puts "    If you intend to send batches with slurm...run..."
puts "        ./run_training.rb --project_name #{base_subdirectory}#{version}_#{timestamp} --fast --slurm"
puts "\n========================================================================\n"
