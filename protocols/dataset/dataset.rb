# TBD:
# verbosity handling
# item selection by name, tag and group
# weighting (name-based, group-based etc)
# spreadsheet-friendly output

require "classes/misc/search_for_file.rb"
require "protocols/dataset/classes/data_set.rb"

module ProtocolDataset
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :batch
	PROTOCOL_DESCRIPTION = "Calculation over a predefined data set and processing of the results"

	INPUT_BLOCKS = [
		InputBlock[:calculation_overwrite, :optional, "Calculation setup that does not fit at the root level"]
	]
	#=======================================================================

	DATASETS_PATH = "data/datasets"
	DATASETS_PATH_DEV = "data/datasets_devel"

	#=======================================================================
	# Classes used for statistical analysis and printing of results
	#=======================================================================
	class DataSetErrors < Array

		def rmse
			sumsq = 0.0
			each{|e| sumsq += e**2 }
			return (sumsq/size)**0.5
		end

		def mue
			sum = 0.0
			each{|e| sum += e.abs }
			return sum/size
		end

		def mse
			sum = 0.0
			each{|e| sum += e }
			return sum/size
		end

		# min, max from array

		def min_abs
			min = 99.9e99
			each{|e|
				min = e.abs if e.abs < min
			}
			return min
		end

		def max_abs
			max = -99.9e99
			each{|e|
				max = e.abs if e.abs > max
			}
			return max
		end
	end

	class DataSetStatistics
		attr_reader :errors
		attr_reader :rel_errors_pct

		def initialize(set_items, jobs, longnums = false, shortnames = false)
			@set_items = set_items
			@jobs = jobs
			
			@statformat = "%10.3f"
			@statformat = "%15.8f" if longnums
			@shortnames = shortnames

			@errors = DataSetErrors.new
			@rel_errors_pct = DataSetErrors.new

			# Average values in the set
			@avg_abs = 0.0
			@avg_ref_abs = 0.0

			@jobs.each_index{|i|
				@errors[i] = @jobs[i].results.energy - @set_items[i].reference_value
				@rel_errors_pct[i] = @errors[i] / @set_items[i].reference_value.abs * 100.0

				@avg_abs += @jobs[i].results.energy.abs
				@avg_ref_abs +=  @set_items[i].reference_value.abs
			}

			@avg_abs /= @jobs.size
			@avg_ref_abs /= @jobs.size
		end

		def print_table
			# header
			puts "=========================================================================================="
			printf("%-40s",'name')
			printf("%10s",'E')
			printf("%10s",'Eref')
			printf("%10s",'error')
			printf("%10s",'error(%)')
			puts
			puts "------------------------------------------------------------------------------------------"

			# Contents
			@set_items.each_index{|i|
				if @shortnames
					name = @set_items[i].shortname
				else
					name = @set_items[i].name
				end
				# Shorten name to 40 characters
				if name.size > 40
					name = name[0..36] + "..."
				end
				printf("%-40s",name)

				printf("%10.3f", @jobs[i].results.energy)
				printf("%10.3f", @set_items[i].reference_value)
				printf("%10.3f", @errors[i])
				printf("%10.3f", @rel_errors_pct[i])
				puts
			}
		end

		def print_stats
			puts "=========================================================================================="
			printf("%-15s#{@statformat}   %s\n", "RMSE", @errors.rmse, "kcal/mol")
			printf("%-15s#{@statformat}   %s\n", "MUE", @errors.mue, "kcal/mol")
			puts "------------------------------------------------------------------------------------------"
			printf("%-15s#{@statformat}   %s\n", "MSE", @errors.mse, "kcal/mol")
			printf("%-15s#{@statformat}   %s\n", "min", @errors.min, "kcal/mol")
			printf("%-15s#{@statformat}   %s\n", "max", @errors.max, "kcal/mol")
			printf("%-15s#{@statformat}   %s\n", "range", @errors.max - @errors.min, "kcal/mol")
			printf("%-15s#{@statformat}   %s\n", "min abs", @errors.min_abs, "kcal/mol")
			printf("%-15s#{@statformat}   %s\n", "max abs", @errors.max_abs, "kcal/mol")
		end

		def print_stats_rel
			puts "=========================================================================================="
			printf("%-15s#{@statformat}   %s\n", "RMSE", @rel_errors_pct.rmse, "%")
			printf("%-15s#{@statformat}   %s\n", "MUE", @rel_errors_pct.mue, "%")
			printf("%-15s#{@statformat}   %s\n", "MSE", @rel_errors_pct.mse, "%")
			printf("%-15s#{@statformat}   %s\n", "min", @rel_errors_pct.min, "%")
			printf("%-15s#{@statformat}   %s\n", "max", @rel_errors_pct.max, "%")
			printf("%-15s#{@statformat}   %s\n", "range", @rel_errors_pct.max - @rel_errors_pct.min, "%")
			printf("%-15s#{@statformat}   %s\n", "min abs", @rel_errors_pct.min_abs, "%")
			printf("%-15s#{@statformat}   %s\n", "max abs", @rel_errors_pct.max_abs, "%")
			puts "------------------------------------------------------------------------------------------"
			printf("%-15s#{@statformat}   %s\n", "MUE/|avg|", @errors.mue / @avg_ref_abs * 100, "%")
		end

		def print_stats_groups(groups)
			puts "=========================================================================================="
			groups.each{|group|
				# Build error set
				gerrors = DataSetErrors.new
				@jobs.each_index{|i|
					if @set_items[i].group && @set_items[i].group.downcase == group.downcase
						gerrors << @jobs[i].results.energy - @set_items[i].reference_value
					end
				}
				# Print
				printf("%-15s%-5s%-5s#{@statformat}   %-5s#{@statformat}   %s\n", group, "(#{gerrors.size})", "RMSE", gerrors.rmse, "MSE", gerrors.mse, "kcal/mol")

			}
		end

		def print_array
			str_array = []
			@set_items.each_index{|i|
				str_array << sprintf("%.3f", @jobs[i].results.energy)
			}
			puts "[ " + str_array.join(", ") + " ]"
		end
	end

	#=======================================================================
	# Protocol methods
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Dataset calculation')
	end

	def prepare_protocol
		# Read dataset
		# Dataset file and directory
		unless fn = SearchForFile.search(:filename => @settings[:dataset], :extensions => ['', '.yaml'], :raise_error => false)
			# The dataset keyword does not point to a file, look for predefined datasets
			dataset_name = @settings[:dataset].strip.gsub(" ","_").downcase
			# Firstly, a yaml file
			# Secondly, file set.yaml in a directory named after the dataset
			unless fn = Cuby.install_dir + "/" + DATASETS_PATH + "/#{dataset_name}.yaml"
				fn = Cuby.install_dir + "/" + DATASETS_PATH + "/#{dataset_name}/set.yaml"
			end

			# Try another location - datsets in development
			unless FileTest.exist?(fn)
				unless fn = Cuby.install_dir + "/" + DATASETS_PATH_DEV + "/#{dataset_name}.yaml"
					fn = Cuby.install_dir + "/" + DATASETS_PATH_DEV + "/#{dataset_name}/set.yaml"
				end
				Cuby::warning("Using a data set under development") if FileTest.exist?(fn)
			end

			Cuby::error("Dataset not found (neither file nor predefined dataset)") unless FileTest.exist?(fn)
		end

		Cuby::error("Dataset keyword does not points to a valid file") unless FileTest.file?(fn)
		@dataset_dir = File.dirname(fn) # Path, later used for loading geometries

		# Load set data from YAML
		File.open(fn, "r") {|f|
			begin
				@set = YAML.load(f)
			rescue
				Cuby::error("Error reading dataset file")
			end
		}
		# Check type
		Cuby::error("Invalid data set file, does not contain DataSet object\n(file #{fn})") unless @set.class == DataSet

		# Sort data set items (used only in preparation of data sets)
		#--------------------------------------------------
		# @set.items.sort!{|a,b|
		# 	a.shortname <=> b.shortname
		# }
		# puts @set.to_yaml
		# exit
		#-------------------------------------------------- 

		# Modify job setup
		if @set.description.dataset_setup
			settings = Settings.from_hash(@set.description.dataset_setup, @settings)
			@settings.copy_from!(settings, :keep)
		end

		# Save original items
		@set.save_orig_items

		# Load items to be modified
		all_items = @set.items

		# Load alternative reference if applicable
		if @settings[:dataset_reference] != "default"
			if @set.alternative_reference
				if @set.alternative_reference.has_key?(@settings[:dataset_reference])
					@set.items.each_index{|i|
						if @set.alternative_reference[@settings[:dataset_reference]][i]
							@set.items[i].reference_value = @set.alternative_reference[@settings[:dataset_reference]][i]
						else
							Cuby::error "Alternative reference '#{@settings[:dataset_reference]}' not found at index #{i}"
						end
					}
				else
					Cuby::error("Alternative reference data requested but not found in the data set file")
				end
			else
				Cuby::error("Alternative reference data requested but not found in the data set file")
			end
		end

		# Add to results data
		if @settings.set?(:dataset_add_ref_to_result)
			if @set.alternative_reference
				if @set.alternative_reference.has_key?(@settings[:dataset_add_ref_to_result])
					@set.items.each_index{|i|
						if @set.alternative_reference[@settings[:dataset_add_ref_to_result]][i]
							@set.items[i].add_to_result = @set.alternative_reference[@settings[:dataset_add_ref_to_result]][i]
						else
							Cuby::error "Alternative reference '#{@settings[:dataset_add_ref_to_result]}' not found at index #{i} (add_ref_to_result)"
						end
					}
				else
					Cuby::error("Alternative reference data requested but not found in the data set file (add_ref_to_result)")
				end
			else
				Cuby::error("Alternative reference data requested but not found in the data set file (add_ref_to_result)")
			end
		end

		# Apply selection filters
		name_filter = Regexp.new(@settings[:dataset_select_name])
		tag_filter = Regexp.new(@settings[:dataset_select_tag])
		name_skip = Regexp.new(@settings[:dataset_skip_name]) if @settings.set?(:dataset_skip_name)
		@set.items = []
		all_items.each{|item|
			if name_filter.match(item.name) && tag_filter.match(item.tags) && (!name_skip || !name_skip.match(item.name))
				@set.items << item
			end
		}

		if @set.items.size == 0
			Cuby::error("There are no item in the set matching the selection criteria")
		end

		# Create child jobs
		@item_jobs = []
		settings_global = Settings.from_hash(@set.description.global_setup, @settings)
		# Check whether it contains job definition
		Cuby::error("Job type for dataset items must be set in description:global_setup in the yaml file") unless settings_global.set?(:job)
		@set.items.each_index{|i|
			item = @set.items[i]
			# Build Settings for the item
			# Copy from root level, exclude the calculation blocks already made and all dataset-specific keywords
			# not needed in the calculations
			settings = @settings.child_block_as_copy(("set_item_" + item.shortname).to_sym, [/^set_item_/, /^dataset.*/, /^calculation_overwrite$/])
			settings.copy_from!(settings_global, :overwrite) # Add global dataset setup from dataset file
			settings.copy_from!(@settings.block(:calculation_overwrite), :overwrite) if @settings.has_block?(:calculation_overwrite) # Add global setup from input
			settings.copy_from!(Settings.from_hash(item.setup, @settings), :overwrite) if item.setup # Add item-specific setup if exist
			# Set the geometry
			if FileTest.exists?(@dataset_dir + "/" + item.geometry)
				settings[:geometry] = @dataset_dir + "/" + item.geometry
			else
				settings[:geometry] = item.geometry
			end

			# Create the jobs
			@item_jobs[i] = Job.new(@name + "_" + item.shortname, settings)
		}

		# Prepare the jobs
		@item_jobs.each{|job| job.prepare }
	end

	def run_protocol(queue)
		# Submit jobs to queue
		@item_jobs.each{|job| job.run}
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			# Get results
			# Modify the results
			if @settings.set?(:dataset_add_ref_to_result)
				@set.items.each_index{|i|
					if @set.items[i].add_to_result
						@item_jobs[i].results.energy += @set.items[i].add_to_result
					else
						raise "This should not happen, check the code!"
					end
				}
			end
			# Processing of results
			case @settings[:dataset_processing]
			when :relative_to_average
				average = 0.0
				@item_jobs.each{|job|
					average += job.results.energy
				}
				average = average / @item_jobs.size
				@item_jobs.each{|job|
					job.results.energy -= average
				}
			end
			# Statistics
			@results = DataSetStatistics.new(@set.items, @item_jobs, @settings[:dataset_long_numbers], @settings[:dataset_short_names])
		end
		return @results
	end

	def print_protocol
		# print dataset info
		# #!#
		#--------------------------------------------------
		# if @set.description.references
		# 	puts "Cite"
		# 	@set.description.references.each_pair{|doi, citation|
		# 		puts "   #{citation} DOI: #{doi}"
		# 	}
		# 	puts
		# end
		#-------------------------------------------------- 
		# Print results
		results.print_table if @settings[:dataset_print].include?(:table)
		results.print_stats if @settings[:dataset_print].include?(:stats)
		results.print_stats_rel if @settings[:dataset_print].include?(:stats_rel)
		results.print_stats_groups(@set.description.groups) if @settings[:dataset_print].include?(:stats_groups) && @set.description.groups
		puts "==========================================================================================" unless @settings[:dataset_print].empty?
		# Save plots
		if @set.plots && @settings[:dataset_save_plots] != :none
			Cuby::log.puts_debug("Saving plots")
			@set.save_plots(@item_jobs, settings)
		end
		# Print array for use in datasets
		results.print_array if @settings[:dataset_print].include?(:results_array)
	end

	def cleanup_protocol
		@item_jobs.each{|job| job.cleanup}
	end

	#=======================================================================
	# Optional protocol functionality
	#=======================================================================
	
	def single_number_result
		return results.errors.rmse
	end
	
end
