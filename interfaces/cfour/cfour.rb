################################################################################
#
# CFOUR interface
#
# Author: Jan Rezac
# Date created: 2013-03-13
# License: Cuby4 license
# Description: Interface for external calculations in CFOUR
# Status: Works
#
################################################################################

#===============================================================================
# Interface to the CFOUR program (mainly coupled clusters calculations)
# Free for academic use (license must be signed)
# http://www.cfour.de/
#===============================================================================

require "classes/tools/output_parser.rb"
require "fileutils"

module InterfaceCfour
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK, more features to be implemented"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :hessian, :ghost_atoms]
	MODIFIER = false
	DIRECTORY = "CFOUR"
	# Methods provided by the interface:
	METHODS = {
		:"hf"		=> [:energy, :gradient, :hessian],
		:"mp2" 		=> [:energy, :gradient, :hessian],
		# :"mp3" 		=> [:energy, :gradient, :hessian],
		# MP3 works but the results differ from other programs
		:"ccsd" 	=> [:energy, :gradient, :hessian],
		:"ccsd(t)" 	=> [:energy, :gradient, :hessian]
	}
	# Unsupported keywords: options common in similar interfaces that can't
	# be used
	UNSUPPORTED_KEYWORDS = [
		:scf_convergence
	]
	#=======================================================================

	# List of Coupled Clusters methods (which share commmon setup)
	CC_METHODS = [:"ccsd", :"ccsd(t)"]

	def prepare_interface
		# Use MRCC
		if @settings[:cfour_force_mrcc]
			@mrcc_used = true
		end

		# Prepare calculation directory
		if calc_dir_mkdir("cfour.in", "cfour.out") == :old
			calc_using_input # Old input should be reused
			return
		end

		# Copy basis set file into it
		cfour_copy_basisset(@settings[:cfour_genbas_file])

		# Write the keyword part of the input
		cfour_write_keywords # This writes file cfour.in

		# Create complete input by combining cfour.in and the current geometry
		cfour_write_input_geo

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		cfour_calculate unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/cfour.out")
		results = cfour_read_results
		return results
	end

	def cleanup_interface
		# Delete large files
		if @settings[:delete_large_files]
			# Call CFOUR's wipe
			system("#{@settings[:cfour_bin_dir]}/xwipeout")
			#!# More files to be deleted manualy
		end

		calc_dir_delete if @settings[:job_cleanup]
	end

	def priority_interface
		int = @geometry.size.to_f
		int *= @settings[:basisset_zeta] if @settings.set?(:basisset_zeta)
		case @settings[:method]
		when :hf
			return int**4
		when :mp2
			return int**5
		when :mp3
			return int**6
		when :"ccsd(t)"
			return int**7
		end
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	# cfour_copy_basisset	
	# cfour_write_keywords
	# cfour_write_input_geo
	# cfour_calculate
	# cfour_read_results

	def cfour_copy_basisset(filename)
		filename = File.expand_path(filename)
		Cuby.log.puts_debug("Cfour interface: using basis set file '#{filename}'")

		unless File.exist?(filename)
			Cuby::error("Cfour basis set file not found at '#{filename}'")
		end

		FileUtils.cp(filename, in_calc_dir("GENBAS"))
	end

	def cfour_write_keywords
		filename = in_calc_dir("cfour.in")
		f = File.open(filename,"w+")
	
		# Begin the calculation setup	
		f.printf "*CFOUR("
		f.puts "COORD=CARTESIAN"
		# Memory requested, in suitable units
		if @settings[:mem] < 10000
			f.puts "MEMSIZE=#{@settings[:mem]}"
			f.puts "MEM_UNIT=MB"
		else
			f.puts "MEMSIZE=#{@settings[:mem] / 1000}"
			f.puts "MEM_UNIT=GB"
		end
		f.puts "BASIS=SPECIAL" # Needed if we want to use ghost atoms
		f.puts "CHARGE=#{@settings[:charge]}"
		f.puts "MULTIPLICITY=#{@settings[:multiplicity]}"
		f.puts "REF=RHF"

		f.puts "ABCDTYPE=AOBASIS" unless @mrcc_used

		# Common setup for coupled clusters
		if CC_METHODS.include?(@settings[:method])
			f.puts "CC_CONV=#{@settings[:correlation_convergence]}"
			# f.puts "CC_EXTRAPOL=NOJACOBI" if @settings[:aces_cc_safe]
			if @mrcc_used
				f.puts "CC_PROGRAM=MRCC"
			else
				f.puts "CC_PROGRAM=ECC"
			end
		end

		case @settings[:method]
		when :hf
			f.puts "CALC=SCF"
		when :mp2
			f.puts "CALC=MP2"
		when :mp3
			f.puts "CALC=MP3"
		when :"ccsd"
			f.puts "CALC=CCSD"
		when :"ccsd(t)"
			f.puts "CALC=CCSD(T)"
		end

		f.puts "FROZEN_CORE=ON" if @settings[:correlation_frozen_core]

		if @settings[:use_symmetry]
			f.puts "SYMMETRY=ON"
		else
			f.puts "SYMMETRY=OFF"
		end

		# Extra keywords from the input
		f.puts @settings[:cfour_keywords] if @settings.set?(:cfour_keywords)

		if what.include?(:hessian)
			f.puts "VIB=EXACT"
		elsif what.include?(:gradient)
			f.puts "DERIVATIVES=FIRST"
		end


		f.puts "SCF_MAXCYC=#{@settings[:scf_cycles]}"

		# Note the ")" at the end, this closes the input section!
		f.puts "SCF_CONV=#{@settings[:density_convergence]})"
		#!# There is no keyword to set the convergence in terms of energy

		f.puts # Blank line

		# Basis set input
		# Using 'special' basis to cover possible use of ghost atoms
		@geometry.each_index{|i|
			f.puts @geometry[i].element.to_s.upcase + ":" + @settings[:basisset].upcase
		}

		f.puts # Blank line

		f.close
	end

	def cfour_write_input_geo
		filename = in_calc_dir("ZMAT")
		f = File.open(filename,"w+")

		# Header
		f.puts "cuby calculation"

		# Geometry
		# Ghost atoms must be named GH
		@geometry.each{|atom|
			element = atom.element.to_s.upcase
			element = 'GH' if atom.properties[:ghost] # Change element for ghost atoms
			f.printf("%-4s %20.12f%20.12f%20.12f\n", element, atom.x, atom.y, atom.z)
		}

		f.puts # Blank line

		# Append the contents of cfour.in
		File.open(in_calc_dir("cfour.in"),"r") {|f_in|
			while s = f_in.gets do
				f.print s
			end
		}

		f.close
	end

	def cfour_calculate
		# Build input with a new geometry
		cfour_write_input_geo

		# Prepare the command to run
		command = "export LD_LIBRARY_PATH=#{@settings[:cfour_ld_path]}:$LD_LIBRARY_PATH;"
		if @settings[:parallel] == 1
			command << "export OMP_NUM_THREADS=1;"
			command << "export PATH=#{@settings[:cfour_bin_dir]}:$PATH;"
		else
			if @settings[:parallel_mode] == :shm
				command << "export OMP_NUM_THREADS=#{@settings[:parallel]};"
				command << "export PATH=#{@settings[:cfour_bin_dir]}:$PATH;"
			else # :mpi or :default
				command << "export PATH=#{@settings[:cfour_bin_dir_mpi]}:$PATH;"
				command << "export CFOUR_NUM_CORES=#{@settings[:parallel]};"
			end
		end
		# MRCC setup
		if @mrcc_used
			command << "export PATH=#{@settings[:cfour_mrcc_bin_dir]}:$PATH;"
		end
		command << "cd #{calc_dir};"
		command << "xwipeout;" # Delete the restart to make gradient calculations possible
		command << "xcfour > cfour.out 2> cfour.err;"

		# Run cfour
		unless system(command)
			Cuby::error "Cfour returned nonzero exit code (calculation #{@name})"
		end
	end

	def cfour_read_hessian_file(filename)
		h = Hessian.zero(@geometry.size * 3)

		f = File.open(in_calc_dir(filename))
		f.gets # Header line

		i = j = 0
		n = @geometry.size * 3
		while line = f.gets
			line.strip.split.map{|x| x.to_f}.each{|x|
				h[i,j] = x * HARTREE2KCAL / BOHR2ANGSTROM / BOHR2ANGSTROM
				i += 1
				if i == n
					i = 0
					j += 1
				end
			}
		end

		f.close
		return h
	end

	def cfour_read_results
		results = Results.new

		parser = OutputParser.new(in_calc_dir("cfour.out"), @name)
		# General output
		parser.add_pattern(:e_scf, /^\s*E\(SCF\)= *(\S+) +/,  {:get => :match_array, :multi => :last})
		# Method-specific output
		case @settings[:method]
		when :mp2
			parser.add_pattern(:corr_mp2, /^\s*E2\(TOT\) *= +(\S+) +/,  {:get => :match_array})
		when :mp3
			parser.add_pattern(:corr_mp2, /^\s*E2\(TOT\) *= +(\S+) +/,  {:get => :match_array})
			parser.add_pattern(:delta_mp3, /^\s*D-MBPT\(3\) +(\S+) +/,  {:get => :match_array})
		when :"ccsd"
			# Gradient:
			parser.add_pattern(:corr_mp2, /^\s*E2\(TOT\) *= +(\S+) +/,  {:get => :match_array})
			if @mrcc_used
				parser.add_pattern(:e_ccsd, /^ Total CCSD energy/)
			else
				parser.add_pattern(:corr_ccsd, /^  CCSD correlation energy/)
			end
		when :"ccsd(t)"
			# Gradient:
			parser.add_pattern(:corr_mp2, /^\s*E2\(TOT\) *= +(\S+) +/,  {:get => :match_array})
			if @mrcc_used
				parser.add_pattern(:e_ccsd, /^ Total CCSD energy/)
				parser.add_pattern(:"e_ccsd(t)", /^ Total CCSD\(T\) energy/)
				parser.add_pattern(:"e_ccsd[t]", /^ Total CCSD\[T\] energy/)
			else
				parser.add_pattern(:e_ccsd, /^  CCSD energy   /)
				parser.add_pattern(:"e_ccsd(t)", /^  CCSD\(T\) energy   /)
				parser.add_pattern(:e_5st, /E5ST to CCSD\(T\)/, {:required => false})
			end
		end

		# Errors
		parser.add_error(/Job has terminated with error flag/, "CFOUR run ended with an error")
		parser.add_error(/Unable to allocate/, "CFOUR complains on insufficient memory")

		# Gradient
		if @what.include?(:gradient)
			grad_new = Gradient.new
			geo_new = Geometry.new

			parser.add_block(:gradient, /Molecular gradient$/, /[0-9]$/, /gradient norm/, {:get => :line_words})
			parser.add_block(:newgeo, /Atomic +Coordinates/, /[0-9]$/, /distance matrix/, {:get => :line_words})
		end

		# The hesian in the output file is already transformed and can not be used
		# for calculations

		# Execute parser
		gi = 0.0
		parser.execute{|name, result, count|
			case name
			when :gradient
				e,n,x,y,z = result
				x = x.to_f * HARTREE2KCAL / BOHR2ANGSTROM
				y = y.to_f * HARTREE2KCAL / BOHR2ANGSTROM
				z = z.to_f * HARTREE2KCAL / BOHR2ANGSTROM
				grad_new << Coordinate[x,y,z]
			when :newgeo
				e,n,x,y,z = result
				x = x.to_f * BOHR2ANGSTROM
				y = y.to_f * BOHR2ANGSTROM
				z = z.to_f * BOHR2ANGSTROM
				geo_new << Atom.new(@geometry[gi].element,x,y,z)
				gi += 1
			end
		}

		if @what.include?(:hessian)
			# Hessian is read from a separate file
			# The one provided in the output file is alraedy projected out and thus incorrect in non-equilibrium geometries
			# The file FCM contains original (symmetrized) Hessian
			# The Hessian must be rotated to the original coordinates
			if @geometry.size == 1
				results.hessian = Hessian.zero(3)
			else
				results.hessian = Hessian.from_calculation_coordinates(@geometry, geo_new, cfour_read_hessian_file("FCM"))
			end
		end

		# Get the total energy and its components
		case @settings[:method]
		when :hf
			results.energy = parser[:e_scf][0] * HARTREE2KCAL
		when :mp2
			results.energy_components[:scf_energy] = parser[:e_scf][0] * HARTREE2KCAL
			results.energy_components[:correlation_mp2] = parser[:corr_mp2][0] * HARTREE2KCAL
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:correlation_mp2]
		when :mp3
			results.energy_components[:scf_energy] = parser[:e_scf][0] * HARTREE2KCAL
			results.energy_components[:correlation_mp2] = parser[:corr_mp2][0] * HARTREE2KCAL
			results.energy_components[:correlation_mp3] = results.energy_components[:correlation_mp2] + parser[:delta_mp3][0] * HARTREE2KCAL
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:correlation_mp3]
		when :"ccsd"
			results.energy_components[:scf_energy] = parser[:e_scf][0] * HARTREE2KCAL
			results.energy_components[:correlation_mp2] = parser[:corr_mp2][0] * HARTREE2KCAL
			if @mrcc_used
				results.energy_components[:"correlation_ccsd"] = parser[:"e_ccsd"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			else
				results.energy_components[:correlation_ccsd] = parser[:corr_ccsd] * HARTREE2KCAL
				results.energy = results.energy_components[:scf_energy] + results.energy_components[:"correlation_ccsd"]
			end
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:"correlation_ccsd"]
		when :"ccsd(t)"
			results.energy_components[:scf_energy] = parser[:e_scf][0] * HARTREE2KCAL
			results.energy_components[:correlation_mp2] = parser[:corr_mp2][0] * HARTREE2KCAL
			if @mrcc_used
				results.energy_components[:"correlation_ccsd"] = parser[:"e_ccsd"] * HARTREE2KCAL - results.energy_components[:scf_energy]
				results.energy_components[:"correlation_ccsd(t)"] = parser[:"e_ccsd(t)"] * HARTREE2KCAL - results.energy_components[:scf_energy]
				results.energy_components[:"correlation_ccsd[t]"] = parser[:"e_ccsd[t]"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			else
				results.energy_components[:correlation_ccsd] = parser[:e_ccsd] * HARTREE2KCAL - results.energy_components[:scf_energy]
				results.energy_components[:"correlation_ccsd(t)"] = parser[:"e_ccsd(t)"] * HARTREE2KCAL - results.energy_components[:scf_energy]
				if parser[:e_5st] # This is nnot printed when gradient is requested
					results.energy_components[:"correlation_ccsd[t]"] = results.energy_components[:"correlation_ccsd(t)"] - parser[:e_5st] * HARTREE2KCAL
				end
			end
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:"correlation_ccsd(t)"]
		end

		# Transform the gradient
		if what.include?(:gradient)
			if @geometry.size == 1
				results.gradient = Gradient.zero(1)
			else
				results.gradient = Gradient.from_calculation_coordinates(@geometry, geo_new, grad_new)
			end
		end


		return results
	end

end
