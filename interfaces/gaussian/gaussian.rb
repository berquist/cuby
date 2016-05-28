################################################################################
#
# Gaussian interface
#
# Author: Jan Rezac
# Date created: 2015-02-02
# License: Cuby4 license
# Description: Interface for external calculations in Gaussian
# Status: Works
#
################################################################################

#===============================================================================
# Interface to the Gaussian package ver. 03 and above
# Commercial software
# http://gaussian.com
#===============================================================================

require "classes/calculation/dft_setup.rb"
require "classes/tools/output_parser.rb"

module InterfaceGaussian
	#=======================================================================
	# Interface header
	#=======================================================================
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "Only simple implementation for now"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy]
	MODIFIER = false
	DIRECTORY = "Gauss"
	# Methods provided by the interface:
	METHODS = {
		:hf		=> [:energy],
		:dft		=> [:energy]
	}
	#SOLVENTS = [:cosmo]
	#ATOMIC_CHARGES = [:mulliken]
	#=======================================================================
	
	INPUT_FILE = "gaussian.com"
	OUTPUT_FILE = "gaussian.log"

	def prepare_interface
		# Things that have to be done even when old input is reused:

		# Prepare calculation directory
		if calc_dir_mkdir("gaussian.header", OUTPUT_FILE) == :old
			calc_using_input # Old input should be reused
			return
		end

		# Write the keyword part of the input
		gaussian_write_keywords # This writes files mopac.header and mopac.footer

		# Create complete input by combining header, footer and the current geometry
		gaussian_write_input_geo

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		gaussian_run unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/"+OUTPUT_FILE)
		return gaussian_read
	end

	def cleanup_interface
		calc_dir_delete if @settings[:job_cleanup]
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def gaussian_write_keywords
		#: Write the header of the input where the calculation is set up
		#: geometry is added later, so that the header can be modified
		#: by the user after the job is prepared

		# HEADER FILE
		filename = in_calc_dir("gaussian.header")
		f = File.open(filename,"w+")

		# Setup
		f.puts "%mem=#{@settings[:mem]}MB"
		f.puts "%Chk=checkpoint.chk" if @settings[:gaussian_checkpoint] != :none

		# Method
		options = ""
		case @settings[:method]
		when :hf
			method = "HF"
		when :dft
			dft_setup = DftSetup.new(interface_dir + "/dft_functionals.yaml", interface_dir + "/dft_grids.yaml")
			method = dft_setup.functional(@settings)
			options += " Integral(Grid=#{dft_setup.grid(@settings)})"
		end

		# Basis set
		method += "/" + @settings[:basisset]
		
		f.puts "#P #{method}#{options} #{@settings[:gaussian_keywords]}"
		f.puts
		f.puts "Gaussian calculation for Cuby"
		f.puts
		f.puts "#{@settings[:charge]} #{@settings[:multiplicity]}"

		f.close

		# FOOTER FILE
		filename = in_calc_dir("gaussian.footer")
		f = File.open(filename,"w+")
		if @settings[:gaussian_extra_block] != ""
			f.puts
			f.puts @settings[:gaussian_extra_block]
			f.puts
		end
		f.puts
		f.puts
		f.close
	end

	def gaussian_write_input_geo
		#: Write final input for mopac, using the previously prepared header and footer
		#: and a current geometry

		# Read header and footer
		header = IO.readlines(in_calc_dir("gaussian.header"))
		footer = IO.readlines(in_calc_dir("gaussian.footer"))

		# Complete input:
		filename = in_calc_dir(INPUT_FILE)
		f = File.open(filename,"w+")

		# Print header from the header file
		f.puts header

		# Write geometry
		@geometry.each_index{|i|
			atom = @geometry[i]
			
			f.printf("%2s  %15.10f %15.10f %15.10f\n",atom.element,atom.x,atom.y,atom.z)
		}

		# Print footer
		f.puts footer

		f.close
	end

	def gaussian_run
		#: Set up environment for gaussian and run it

		# Build input with the current geometry
		gaussian_write_input_geo

		# Run gaussian
		command =  "export GAUSS_EXEDIR=\"#{@settings[:gaussian_bin_dir]}\";"
		command << "export LD_LIBRARY_PATH=\"#{@settings[:gaussian_bin_dir]}\":$LD_LIBRARY_PATH;"
		command << "export GAUSS_SCRDIR=\"/tmp/gaussian/gaussian_$$\";"
		command << "cd #{calc_dir};"
		command << "mkdir -p $GAUSS_SCRDIR;"
		command << "$GAUSS_EXEDIR/#{@settings[:gaussian_version]} #{INPUT_FILE} 2> gaussian.err > gaussian.stdout;"
		command << "rm -rf $GAUSS_SCRDIR;"
		unless system(command)
			Cuby::error "Gaussian returned nonzero exit code (calculation #{@name})"
		end

		# Error check
		Cuby::error "Gaussian output file not found" unless FileTest.exist?(in_calc_dir(OUTPUT_FILE))

		# Error check
		f = File.open(in_calc_dir(OUTPUT_FILE)) 	# Get last line of the file:
		f.seek(-100, IO::SEEK_END)			# Set position 100 bytes before EOF
		f.gets						
		lastline = f.gets
		f.close
		Cuby::error "Gaussian run did not terminate normally" unless lastline =~ /Normal termination of Gaussian/

		# Convert checkpoint to text
		if @settings[:gaussian_checkpoint] == :text
			system("cd #{calc_dir}; #{@settings[:gaussian_bin_dir]}/formchk checkpoint.chk 2>&1 > formchk.out")
		end
	end

	def gaussian_read # => Results
		# Initialize results
		results = Results.new

		# Try reading the results block first
		results_block = gaussian_get_results_block

		# Debugging: print the results block
		if Cuby::log.logs[0].verbosity == :debug
			Cuby::log.puts_debug "Gaussian results archive:"
			results_block.each_index{|i|
				Cuby::log.puts_debug "block #{i}:" 
				Cuby::log.puts_debug results_block[i]
			}
			Cuby::log.puts_debug ""
		end

		case @settings[:method]
		when :hf
			# Energy from the results block
			results.energy = results_block[4][2].gsub("HF=","").to_f * HARTREE2KCAL
		when :dft
			# Energy from the results block
			results.energy = results_block[4][2].gsub("HF=","").to_f * HARTREE2KCAL
		end


		return results
	end

	def gaussian_get_results_block
		#: Find the machine-readable "archive" block and return it as an array
		f = File.open(in_calc_dir(OUTPUT_FILE), "r")
		results_s = ""
		save_this = false
		f.each_line {|line|
			if line =~ /^ *1\\1/
				save_this = true
			end
			results_s << line.chomp.strip if save_this
			if line =~ /\\\\@$/
				save_this = false
			end
		}
		f.close
		a = results_s.split("\\\\").collect{|b| b.split("\\")}

		return a
	end

end
