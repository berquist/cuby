################################################################################
#
# MOLPRO interface
#
# Author: Jan Rezac
# Date created: 2014-01-20
# License: Cuby4 license
# Description: Interface for external calculations using MOLPRO
# Status: Works
#
################################################################################

#===============================================================================
# Interface to MOLPRO
# http://www.molpro.net/
#===============================================================================

# TBD (compared to cuby3)
# * mp3
# * dft
# * mp2c
# * mp2-f12, ccsd-f12, ccsd(t)-f12
# * lmp2
# * sapt
# * basis set from file
# * basisset_elements
# * auxiliary basis choice
# * point charges
# * open shell calculations

require "erb" # ERB templating system, used for construction of the input
require "classes/tools/output_parser.rb"

module InterfaceMolpro
	#=======================================================================
	# Interface header
	#=======================================================================
	DEVELOPMENT_FLAG = :warning
	DEVELOPMENT_STATUS = "Interface works but was not tested yet"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :ghost_atoms, :point_charges]
	MODIFIER = false
	DIRECTORY = "MOLPRO"
	# Methods provided by the interface:
	METHODS = {
		:"hf"		=> [:energy, :gradient, :point_charges],
		:"mp2"		=> [:energy, :gradient, :point_charges],
		:"ccsd"		=> [:energy],
		:"ccsd(t)"	=> [:energy],
		:"mp3"		=> [:energy],
		#:"dft"		=> [:energy, :gradient],
	}
	SOLVENTS = [:cosmo]
	#=======================================================================
	
	MOLPRO_INPUT =  "molpro.com"
	MOLPRO_OUTPUT = "molpro.out"

	def prepare_interface
		# Check setup
		molpro_check_setup

		# Prepare calculation directory
		if calc_dir_mkdir(MOLPRO_INPUT, MOLPRO_OUTPUT) == :old
			calc_using_input # Old input should be reused
			return
		end

		# Write input
		molpro_write_input

		# Write geometry
		@geometry.write_xyz(:file => in_calc_dir("geo.xyz"))

		# Write point charges
		molpro_write_charges if @point_charges

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		# Write current geometry
		@geometry.write_xyz(:file => in_calc_dir("geo.xyz"))

		# Write point charges
		molpro_write_charges if @point_charges

		# Run molpro
		molpro_run unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(in_calc_dir(MOLPRO_OUTPUT))

		# Read results
		return molpro_read
	end

	def cleanup_interface
		calc_dir_delete if @settings[:job_cleanup]
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def molpro_check_setup
		if @settings[:density_fitting] == :scf
			Cuby::error "Correlation energy calculation without density fittingcan not be performed\nafter a SCF calculation with density fitting"
		end

		if @settings[:method] == :mp2 && [:both, :correlation].include?(@settings[:density_fitting]) && @what.include?(:gradient)
			Cuby::error("Gradient not available for MP2 with density fitting in Molpro")
		end

		if @settings[:method] == :mp3 && @settings[:density_fitting] != :none
			Cuby::error("MP3 with density fitting not available in Molpro")
		end

		if @settings[:spin_restricted] == false && @settings[:density_fitting] != :none
			Cuby::error("Spin-unrestricted calculations not available with desnity fitting")
		end

		if @settings[:spin_restricted] == false && ! [:hf, :mp2].include?(@settings[:method])
			Cuby::error("Spin-unrestricted calculations available only for HF and MP2")
		end

		if @settings[:method] == :mp3 && @settings[:multiplicity] != 1
			Cuby::error("Open-shell MP3 calculations not available")
		end
	end
	
	def molpro_write_input
		# Generate the input file from the templates
		file = File.open(in_calc_dir(MOLPRO_INPUT),"w+")

		# Open shell
		if @settings[:spin_restricted]
			if @settings[:multiplicity] != 1
				openshell = "r"
			else
				openshell = ""
			end
		else
			openshell = "u"
		end

		# General setup
		template = ERB.new(IO.read(interface_dir + "/templates/general.erb"))
		file.print(template.result(binding))

		# HF setup
		template = ERB.new(IO.read(interface_dir + "/templates/hf.erb"))
		file.print(template.result(binding))

		# Method-specific templates
		case @settings[:method]
		when :mp2
			# MP2 setup
			template = ERB.new(IO.read(interface_dir + "/templates/mp2.erb"))
			file.print(template.result(binding))
		when :mp3
			# MP3 setup
			template = ERB.new(IO.read(interface_dir + "/templates/mp3.erb"))
			file.print(template.result(binding))
		when :"ccsd"
			# CCSD setup
			template = ERB.new(IO.read(interface_dir + "/templates/ccsd.erb"))
			file.print(template.result(binding))
		when :"ccsd(t)"
			# CCSD(T) setup
			template = ERB.new(IO.read(interface_dir + "/templates/ccsdt.erb"))
			file.print(template.result(binding))
		end

		file.close
	end

	def molpro_write_charges
		f = File.open(in_calc_dir("extcharges.txt"), "w+")
		f.puts "Charges"
		f.puts @point_charges.size
		@point_charges.each{|pch|
			f.puts "#{pch.x},#{pch.y},#{pch.z},#{pch.charge},1"
		}
		f.close
	end

	def molpro_run

		# delete molpro.out and molpro.xml
		["molpro.out", "molpro.xml"].each{|delfile|
			File.delete("#{calc_dir}/#{delfile}") if FileTest.exist?("#{calc_dir}/#{delfile}")
		}

		# Run calculation
		command = "cd #{calc_dir};"
		command += "export OMP_NUM_THREADS=#{@settings[:parallel]};"
		command += "#{@settings[:molpro_exe]} -n #{@settings[:parallel]}"
	       	command += " -L #{@settings[:molpro_lib_dir]}" if @settings.set?(:molpro_lib_dir)
		command += " #{MOLPRO_INPUT} > molpro.err 2>&1"
		unless system(command)
			# Look for common errors in the output
			parser = OutputParser.new(in_calc_dir(MOLPRO_OUTPUT), @name)
			parser.add_error(/Cannot find default basis/,				"MOLPRO failed, unknown default basis set")
			parser.add_error(/Basis library exhausted/,				"MOLPRO failed, basis set not found")
			parser.add_error(/Unknown parameter .* for set HFOPT:DFUNC:KSOPT/,	"MOLPRO failed, unknown DFT functional")
			parser.add_error(/SPIN not possible/,					"MOLPRO failed, wrong multiplicity/charge combination")
			parser.add_error(/NOT ENOUGH MEMORY/,					"MOLPRO failed, not enough memory")
			parser.add_error(/words of memory are needed/,				"MOLPRO failed, not enough memory")
			parser.add_error(/Catastrophic lack of memory/,				"MOLPRO failed, not enough memory")
			parser.add_error(/insufficient memory available/,			"MOLPRO failed, not enough memory")
			parser.add_error(/I\/O error/,						"MOLPRO failed, I/O error")
			parser.execute
			# Nonzero exit code:
			# no erro for now, it seems that molpro returns it even when everything is OK
		end
	end

	def molpro_read
		results = Results.new
		parser = OutputParser.new(in_calc_dir(MOLPRO_OUTPUT), @name)

		# Open shell
		if @settings[:spin_restricted]
			if @settings[:multiplicity] != 1
				openshell = "R"
			else
				openshell = ""
			end
		else
			openshell = "U"
		end

		# SCF energy
		parser.add_pattern(:e_scf, /^\s*SETTING ENERGYSCF += +([^ ]*)/, :get => :submatch)

		# Correlation energy
		if @settings[:method] == :mp2
			if @settings[:spin_restricted]
				parser.add_pattern(:corr_mp2,  /^\s*[DF-]*#{openshell}MP2 correlation energy\s+([^ ]*)/, :get => :submatch)
				if @settings[:multiplicity] == 1
					parser.add_pattern(:corr_mp2s, /^\s*MP2 singlet pair energy\s+([^ ]*)/, :get => :submatch)
					parser.add_pattern(:corr_mp2t, /^\s*MP2 triplet pair energy\s+([^ ]*)/, :get => :submatch)
				end
			else
				parser.add_pattern(:corr_mp2,  /^\s*Correlation energy\s+([^ ]*)/, :get => :submatch)
			end
		end

		if @settings[:method] == :mp3
			parser.add_pattern(:corr_mp2,  /^\s*MP2:\s+([^ ]*)/, :get => :submatch)
			parser.add_pattern(:corr_mp3,  /^\s*MP3 correlation energy\s+([^ ]*)/, :get => :submatch)
			parser.add_pattern(:corr_mp3s, /^\s*MP3 singlet pair energy\s+([^ ]*)/, :get => :submatch)
			parser.add_pattern(:corr_mp3t, /^\s*MP3 triplet pair energy\s+([^ ]*)/, :get => :submatch)
		end

		if [:ccsd, :"ccsd(t)"].include? @settings[:method]
			if @settings[:multiplicity] == 1
				# MP2
				parser.add_pattern(:corr_mp2,  /^\s*MP2 correlation energy:\s+([^ ]*)/, :get => :submatch)
				parser.add_pattern(:corr_mp2s, /^\s*MP2 singlet pair energy:\s+([^ ]*)/, :get => :submatch)
				parser.add_pattern(:corr_mp2t, /^\s*MP2 triplet pair energy:\s+([^ ]*)/, :get => :submatch)
				# CCSD
				parser.add_pattern(:corr_ccsd,  /^\s*CCSD correlation energy\s+([^ ]*)/, :get => :submatch)
				parser.add_pattern(:corr_ccsd_s,  /^\s*CCSD singlet pair energy\s+([^ ]*)/, :get => :submatch)
				parser.add_pattern(:corr_ccsd_t,  /^\s*CCSD triplet pair energy\s+([^ ]*)/, :get => :submatch)
			else
				# MP2
				parser.add_pattern(:corr_mp2,  /^\s*RHF-RMP2 correlation energy\s+([^ ]*)/, :get => :submatch)
				# CCSD
				parser.add_pattern(:corr_ccsd,  /^\s*RCCSD correlation energy\s+([^ ]*)/, :get => :submatch)
			end
		end

		if @settings[:method] == :"ccsd(t)"
			# (T) and [T]
			parser.add_pattern(:corr_triples,  /^\s*Triples \(T\) contribution\s+([^ ]*)/, :get => :submatch)
			parser.add_pattern(:"energy_ccsd[t]",  /^\s*.*CCSD\[T\] energy\s+([^ ]*)/, :get => :submatch)
		end

		# Gradient
		if @what.include?(:gradient)
			case @settings[:method]
			when :hf
				grad_type = "SCF"
			when :mp2
				grad_type = "MP2"
			end
			results.gradient = Gradient.new
			parser.add_block(
				:gradient,
				/#{grad_type} GRADIENT FOR STATE/, 
				/^ +[0-9]+ +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+)/, 
				/Nuclear/,
				{:get => :match_array, :type => :float}
			)
		end

		parser.execute{|name, result, count|
			case name
			when :gradient
				x,y,z = result
				results.gradient << Coordinate[x,y,z] * (HARTREE2KCAL / BOHR2ANGSTROM)
			end
		}

		# SCF energy
		results.energy = parser[:e_scf] * HARTREE2KCAL
		results.energy_components[:scf_energy] = results.energy

		# Evaluate MP2 components
		if [:mp2, :ccsd, :"ccsd(t)"].include? @settings[:method]
			results.energy_components[:correlation_mp2] = parser[:corr_mp2] * HARTREE2KCAL
			# Spin components only for singlet
			if @settings[:multiplicity] == 1 && @settings[:spin_restricted]
				# Correct same-spin and other-spin components have to be calculated
				ss = parser[:corr_mp2t] * HARTREE2KCAL / 1.5
				results.energy_components[:correlation_mp2_os] = results.energy_components[:correlation_mp2] - ss
				results.energy_components[:correlation_mp2_ss] = ss
			end
		end

		if @settings[:method] == :mp3
			results.energy_components[:correlation_mp2] = parser[:corr_mp2] * HARTREE2KCAL
			results.energy_components[:correlation_mp3] = parser[:corr_mp3] * HARTREE2KCAL
			# Correct same-spin and other-spin components have to be calculated
			ss = parser[:corr_mp3t] * HARTREE2KCAL / 1.5
			results.energy_components[:correlation_mp3_os] = results.energy_components[:correlation_mp3] - ss
			results.energy_components[:correlation_mp3_ss] = ss
		end

		# CCSD components
		if [:ccsd, :"ccsd(t)"].include? @settings[:method]
			results.energy_components[:correlation_ccsd] = parser[:corr_ccsd] * HARTREE2KCAL
			# Correct same-spin and other-spin components have to be calculated
			if @settings[:multiplicity] == 1 && @settings[:spin_restricted]
				ss = parser[:corr_ccsd_t] * HARTREE2KCAL / 1.5
				results.energy_components[:correlation_ccsd_os] = results.energy_components[:correlation_ccsd] - ss
				results.energy_components[:correlation_ccsd_ss] = ss
			end
		end

		# Post-HF
		case @settings[:method]
		when :mp2
			# Add correlation to the SCF energy
			results.energy += results.energy_components[:correlation_mp2]
		when :ccsd
			# Add correlation to the SCF energy
			results.energy += results.energy_components[:correlation_ccsd]
		when :"ccsd(t)"
			results.energy_components[:"correlation_ccsd[t]"] = parser[:"energy_ccsd[t]"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			results.energy_components[:"correlation_ccsd(t)"] = results.energy_components[:correlation_ccsd] + parser[:corr_triples] * HARTREE2KCAL
			# Add correlation to the SCF energy
			results.energy += results.energy_components[:"correlation_ccsd(t)"]
		end

		return results
	end
end
