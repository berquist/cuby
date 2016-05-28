################################################################################
#
# MRCC interface
#
# Author: Jan Rezac
# Date created: 2015-06-05
# License: Cuby4 license
# Description: Interface for external calculations in MRCC
# Status: Works
#
################################################################################

#===============================================================================
# Interface to the MRCC program (mainly coupled clusters calculations)
# Free for academic use (license must be signed)
# http://www.mrcc.hu
#===============================================================================

require "classes/tools/output_parser.rb"
require "fileutils"

module InterfaceMrcc
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK, more features to be implemented"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :ghost_atoms]
	MODIFIER = false
	DIRECTORY = "MRCC"
	# Methods provided by the interface:
	METHODS = {
		:"hf"		=> [:energy],
		:"mp2" 		=> [:energy],
		:"ccsd" 	=> [:energy],
		:"ccsd(t)" 	=> [:energy]
	}
	#=======================================================================

	# List of Coupled Clusters methods (which share commmon setup)
	CC_METHODS = [:"ccsd", :"ccsd(t)"]

	def prepare_interface
		# Prepare calculation directory
		if calc_dir_mkdir("MINP", "mrcc.out") == :old
			calc_using_input # Old input should be reused
			return
		end

		# Write the keyword part of the input
		mrcc_write_keywords

		# Create complete input by combining mrcc.in and the current geometry
		mrcc_write_input_geo

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		mrcc_calculate unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/mrcc.out")
		if @settings[:method] == :hf || @settings[:method] == :mp2
			results = mrcc_read_results_iface
		else
			results = mrcc_read_results_output
		end

		return results
	end

	def cleanup_interface
		# Delete large files
		if @settings[:job_cleanup]
			calc_dir_delete 
		elsif @settings[:delete_large_files]
			# Delete all but input and important output
			system("mv #{calc_dir} #{calc_dir}_2; mkdir #{calc_dir}; cd #{calc_dir}_2; mv mrcc.* iface _calc_setup COORD.xyz ../#{calc_dir}; cd ..; rm -rf #{calc_dir}_2")
		end
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

	MRCC_METHODS = {
		:hf => "HF",
		:mp2 => "MP2",
		:ccsd => "CCSD",
		:"ccsd(t)" => "CCSD(T)"
	}
	def mrcc_write_keywords
		filename = in_calc_dir("mrcc.in")
		f = File.open(filename,"w+")

		#f.puts "=#{@settings[:]}"
	
		# General setup
		f.puts "# MRCC calculation for Cuby"
		f.puts "# General setup"
		f.puts "mem=#{@settings[:mem]}MB"

		# Density fitting
		df = 'DF-'
		unless @settings[:density_fitting] == :both
			Cuby::warning "Density fitting is always on in MRCC calculations"
		end

		# Method
		f.puts "# Method setup"
		f.puts "calc=" + df + MRCC_METHODS[@settings[:method]]

		# Local CC
		if @settings[:mrcc_local_cc] == :cim
			# Cluster-in-molecules
			# As described in http://dx.doi.org/10.1063/1.3632085
			f.puts "localcc=on"
			f.puts "lmp2dens=on"
			f.puts "domrad=inf"
			#lnoepso
			#lnoepsv
			#
			raise "LCC/CIM setup not finished"
		elsif @settings[:mrcc_local_cc] == :lno
			# Local natural orbitals
			# Setup from http://dx.doi.org/10.1063/1.4819401
			# Confirmed by M. Kallay
			f.puts "localcc=on"	# Switches Local CC on
			f.puts "lmp2dens=on" 	# Default value
			f.puts "domrad=#{@settings[:mrcc_lno_domrad]}"	# Default value 10
			f.puts "bpcompo=0.985"	# Default value, in the original paper 0.98 was used
			f.puts "lnoepso=1e-#{@settings[:mrcc_lno_eps_o]}"	# Threshold determining accuracy, default 5
			f.puts "lnoepsv=1e-#{@settings[:mrcc_lno_eps_v]}"	# Threshold determining accuracy, default 6
		end

		# Thresholds
		f.puts "# Thresholds"
		f.puts "scftol=#{@settings[:scf_convergence]}"
		f.puts "cctol=#{@settings[:correlation_convergence]}"

		# Basisset
		f.puts "# Basis set"
		f.puts "basis="+@settings[:basisset]

		# System setup
		f.puts "# System setup"
		f.puts "charge=#{@settings[:charge]}"
		f.puts @settings[:correlation_frozen_core] ? "core=frozen" : "core=0"

		# Ghost atoms
		ghosts = []
		@geometry.each_index{|i|
			ghosts << (i + 1) if @geometry[i].properties[:ghost]
		}
		if ghosts.size > 0
			f.puts "# Ghost atoms"
			f.puts "ghost=serialno"
			f.puts ghosts.join(",")
		end
			

		# Geometry header
		f.puts "# Geometry"
		f.puts "unit=angs"
		f.puts "geom=xyz"

		f.close
	end

	def mrcc_write_input_geo
		filename = in_calc_dir("MINP")
		f = File.open(filename,"w+")

		# Insert the contents of mrcc.in
		File.open(in_calc_dir("mrcc.in"),"r") {|f_in|
			while s = f_in.gets do
				f.print s
			end
		}

		# Add geometry
		@geometry.write_xyz(:file => f)

		f.close
	end

	def mrcc_calculate
		# Build input with a new geometry
		mrcc_write_input_geo

		# Prepare commandline to run
		command = "cd #{calc_dir};"
		command << "export PATH=#{@settings[:mrcc_bin_dir]}:$PATH;"
		command << "export OMP_NUM_THREADS=#{@settings[:parallel]};"
		command << "dmrcc > mrcc.out 2> mrcc.err;"

		# Run mrcc
		unless system(command)
			Cuby::error "MRCC returned nonzero exit code (calculation #{@name})"
		end
	end

	def mrcc_read_results_iface
		results = Results.new

		parser = OutputParser.new(in_calc_dir("iface"), @name)

		# SCF
		parser.add_pattern(:e_scf, /^ENERGY *SCF *\S+ *\S+ *\S+  ([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => true})
		# MP2
		parser.add_pattern(:e_mp2, /^ENERGY *MP2 +\S+ *\S+ *\S+  ([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => true})
		parser.add_pattern(:corr_mp2s, /^ENERGY *MP2singlet +\S+ *\S+ *\S+  ([-]*[0-9]*\.[0-9,E,\+,-]+)/,   {:get => :submatch, :type => :float, :required => true})
		parser.add_pattern(:corr_mp2t, /^ENERGY *MP2triplet +\S+ *\S+ *\S+  ([-]*[0-9]*\.[0-9,E,\+,-]+)/,   {:get => :submatch, :type => :float, :required => true})

		# Execute parser
		parser.execute


		# Get the total energy and its components
		case @settings[:method]
		when :hf
			results.energy = parser[:e_scf] * HARTREE2KCAL
		when :mp2
			results.energy_components[:scf_energy] = parser[:e_scf] * HARTREE2KCAL
			results.energy = parser[:e_mp2] * HARTREE2KCAL
			results.energy_components[:correlation_mp2] = results.energy - results.energy_components[:scf_energy]
			# Spin components:
			results.energy_components[:correlation_mp2_os] = parser[:corr_mp2s] * HARTREE2KCAL
			results.energy_components[:correlation_mp2_ss] = parser[:corr_mp2t] * HARTREE2KCAL
		end

		return results
	end

	def mrcc_read_results_output
		results = Results.new

		parser = OutputParser.new(in_calc_dir("mrcc.out"), @name)

		# SCF
		parser.add_pattern(:e_scf, /^ *\*\*\*FINAL HARTREE-FOCK ENERGY: +([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => true})

		if @settings[:mrcc_local_cc] == :none
			# CCSD
			parser.add_pattern(:e_mp2, /^ *MP2 energy \[au\]: +([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => false})
			parser.add_pattern(:corr_ccsd, /^ *CCSD correlation energy \[au\]: +([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => false})
			# CCSD(T)
			parser.add_pattern(:"corr_ccsd(t)", /^ *CCSD\(T\) correlation energy \[au\]: +([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => false})
		else
			# Local CC - read only final local values
			# MP2 correction
			parser.add_pattern(:local_correction, /^ *Local MP2 correction \[au\]: +([-]*[0-9]*\.[0-9,E,\+]+)/, {:get => :submatch, :type => :float})
			# CCSD
			parser.add_pattern(:e_mp2, /^ *Total LMP2 energy \[au\]: +([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => false})
			parser.add_pattern(:corr_ccsd, /^ *CCSD correlation energy \[au\]: +([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => false, :multi => :last})
			# CCSD(T)
			parser.add_pattern(:"corr_ccsd(t)", /^ *CCSD\(T\) correlation energy \[au\]: +([-]*[0-9]*\.[0-9,E,\+]+)/,   {:get => :submatch, :type => :float, :required => false, :multi => :last})
		end

		# Execute parser
		parser.execute

		# Get the total energy and its components
		if [:ccsd, :"ccsd(t)"].include?(@settings[:method])
			results.energy_components[:scf_energy] = parser[:e_scf] * HARTREE2KCAL
			results.energy_components[:correlation_mp2] = parser[:e_mp2] * HARTREE2KCAL - results.energy_components[:scf_energy]
			results.energy_components[:correlation_ccsd] = parser[:corr_ccsd] * HARTREE2KCAL
		end
		if @settings[:method] == :ccsd
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:correlation_ccsd]
		end
		if @settings[:method] == :"ccsd(t)"
			results.energy_components[:"correlation_ccsd(t)"] = parser[:"corr_ccsd(t)"] * HARTREE2KCAL
			results.energy = results.energy_components[:scf_energy] +  results.energy_components[:"correlation_ccsd(t)"]
		end

		# Local CC
		if @settings[:mrcc_local_cc] != :none
			results.energy_components[:local_mp2_correction] = parser[:local_correction] * HARTREE2KCAL
		end

		return results
	end

end
