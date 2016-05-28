################################################################################
#
# Turbomole interface
#
# Author: Jan Rezac
# Date created: 2012-09-25
# License: Cuby4 license
# Description: Interface for external calculations in Turbomole
# Status: Works
#
################################################################################

#===============================================================================
# Interface to the Turbomole program
# Commercial
# http://www.turbomole-gmbh.com/
#===============================================================================

#===============================================================================
# ToDo
#-------------------------------------------------------------------------------
# ( ) Auxiliary basis set fix when uknown basis set is used
# ( ) Reuse previous MOs
# ( ) Save energy components
# ( ) MP2 without RI
# ( ) CCSD
# ( ) Basis sets from basis set exchange
#===============================================================================


require "erb" # ERB templating system, used for construction of the input
require "classes/tools/output_parser.rb"
require "classes/misc/search_for_file.rb"
require "classes/misc/grep.rb"

module InterfaceTurbomole
	#=======================================================================
	# Interface header
	#=======================================================================
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK, more features to be implemented"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :ghost_atoms, :electric_field, :hessian, :atomic_charges, :point_charges, :point_charges_gradient, :mos]
	MODIFIER = false
	DIRECTORY = "TMOLE"
	# Methods provided by the interface:
	METHODS = {
		:"hf"		=> [:energy, :gradient, :hessian, :atomic_charges, :point_charges, :point_charges_gradient, :mos],
		:"mp2" 		=> [:energy, :gradient],
		:"mp3" 		=> [:energy],
		:"dft" 		=> [:energy, :gradient, :hessian, :atomic_charges, :point_charges, :point_charges_gradient, :mos],
		:"ccsd" 	=> [:energy],
		:"ccsd(t)" 	=> [:energy],
		:"rpa" 		=> [:energy],
	}
	SOLVENTS = [:cosmo]
	# Atomic charges
	ATOMIC_CHARGES = [:nbo, :mulliken]
	#=======================================================================
	
	# Extra output files
	FILE_POLARIZABILITY = "polarizability"


	def prepare_interface
		# Check version and compatibility

		# Multiplicity
		if @settings[:multiplicity] != 1 && @settings[:spin_restricted]
			Cuby::error "For spin multiplicity > 1, only unrestricted calculations are supported,\n (set keyword 'spin_restricted: no')"
		end
		
		# RI
		@scf_ri = @settings[:density_fitting] == :scf || @settings[:density_fitting] == :both
		@correlation_ri = @settings[:density_fitting] == :correlation || @settings[:density_fitting] == :both

		if @scf_ri && @what.include?(:hessian)
			Cuby::error("Hessian calculation can not be performed with the RI approximation")
		end

		# Method checks
		if @settings[:method] == :mp2 && !@correlation_ri
			# MP2 not implemented, only RI-MP2 works
			Cuby::error("MP2 calculation without RI currently not supported by the interface", self)
		end
		if @settings[:method] == :rpa && (!@correlation_ri || !@scf_ri)
			# Only RI-RPA possible
			Cuby::error("RPA calculation without RI not supported", self)
		end
		if @settings[:spin_component_scaling] && ! (@settings[:method] == :mp2 && @correlation_ri)
			# SCS available only for RI-MP2
			Cuby::error("SCS is at present available only for RI-MP2")
		end
		if @settings[:method] == :mp3 && @settings[:turbomole_version] < 6.5
			Cuby::error("MP3 is not available in turbomole version lower than 6.5")
		end

		if @settings[:explicit_correlation] != :none
			unless @settings[:method] == :mp2 && ! @what.include?(:gradient)
				Cuby::error("Explicit corrlation is available only for MP2 energy calculations for now")
			end
		end

		# Symmetry
		if @settings[:use_symmetry]
			 if @settings[:explicit_correlation] != :none
				Cuby::error "F12 calculations possible only without use_symmetry (#{@name})"
			end
			 @geometry.each{|atom|
				if atom.properties[:ghost]
					Cuby::warning "Symmetry can't be used with ghost atoms in turbomole,\nswitching it off (#{@name})"
					@settings[:use_symmetry] = false
				end
			 }
		end

		# Prepare calculation directory
		if calc_dir_mkdir("control", "energy") == :old
			calc_using_input
			return
		end

		# Write geometry
		GeoFile::TurbomoleCoord.write(@geometry, :file => calc_dir+"/coord")

		# Prepare input
		calc_writing_input
		turbomole_run_define
		turbomole_cosmoprep if @settings[:solvent_model] == :cosmo
		turbomole_fix_control
	end

	def calculate_interface
		turbomole_run_calculation unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/energy")

		results = turbomole_read_results

		return results
	end

	def cleanup_interface
		if @settings[:job_cleanup]
			# Delete everything
			calc_dir_delete 
		elsif @settings[:delete_large_files]
			# Delete only large files (and some unnecessary files as well)
			system("cd #{calc_dir}; rm -f mos CC* slave*.output master ewald_1 statistics")
		end
	end

	def priority_interface
		int = @geometry.size.to_f
		int *= @settings[:basisset_zeta] if @settings.set?(:basisset_zeta)
		case @settings[:method]
		when :dft
			return int**4
		when :hf
			return int**4
		when :mp2
			return int**5
		when :mp3
			return int**6
		when :"ccsd"
			return int**6
		when :"ccsd(t)"
			return int**7
		end
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def turbomole_run_define

		# Copy basis set file if needed
		if @settings.set?(:basisset_file)
			bfn = SearchForFile.search(:filename => @settings[:basisset_file])
			FileUtils.cp(bfn, calc_dir+"/basis")
		end

		# Load input templates: general
		template_common = ERB.new(IO.read(interface_dir + "/templates/define_common.erb"))

		# Create input for define
		File.open(calc_dir+"/define.in","w+") {|file|
			file.print(template_common.result(binding))
		}

		# Copy the control_template file
		FileUtils.cp(interface_dir+"/templates/control_template", calc_dir+"/control_template")

		# Run define to do the general setup
		system("cd #{calc_dir}; export TURBODIR=#{@settings[:turbomole_turbodir]}; cat define.in | #{@settings[:turbomole_bin_dir]}/define > define.out 2>&1");

		# Look for errors in define output
		unless Grep::file_unicode(calc_dir+"/define.out", /define ended normally/)
			Cuby::error("Turbomole define script failed, see #{calc_dir}/define.out for details")
		end

		# Check for specific errors / define output
		parser = OutputParser.new(in_calc_dir("define.out"), @name)
		parser.add_error(/THERE ARE NO DATA SETS CATALOGUED IN FILE/, "Basis set not found in the library")
		parser.add_pattern(:multiplicity, /^ *FOUND .* WITH MULTIPLICITY +([0-9]+)/, :get => :submatch, :type => :as_is, :required => false)
		parser.add_pattern(:multiplicity_cs, /^ *FOUND CLOSED SHELL SYSTEM/, :get => :line, :type => :as_is, :required => false)
		parser.add_pattern(:electrons_in_atom, /^ *NUMBER OF ELECTRONS IN YOUR ATOM IS +([0-9]+)/, :get => :submatch, :type => :as_is, :required => false)
		parser.add_pattern(:multiplicity_s_trivial, /SINCE ALL MOLECULAR ORBITALS ARE DOUBLY/, :get => :line, :type => :as_is, :required => false)
		parser.execute

		# Check multiplicity
		if @settings[:multiplicity] == 1
			multiplicity = -1
			if parser[:multiplicity_cs] || parser[:multiplicity_s_trivial]
				multiplicity = 1 
			elsif parser[:multiplicity]
				multiplicity = parser[:multiplicity].to_i
			elsif parser[:electrons_in_atom]
				# It is an atom, multiplicity not in the output of define
				# -> not checked, assume the input is correct
				multiplicity = @settings[:multiplicity]
			end
			if multiplicity != @settings[:multiplicity]
				Cuby::error("Multiplicity determined by turbomole's define is different from the one\nspecified in the input. (#{@name})")
			end
		end


		# Memory setup, depends on parallelization mode
		# The :mem keyword should specify total memory, per-core memory is calculated here
		# ri_scf_core and ri_corr_core variables have to be set

		# In serial jobs, MPI and GA parallelization, memory per core is used
		if @settings.set?(:mem_core_scf)
			ri_scf_core =  @settings[:mem] * @settings[:mem_core_scf] / 100 / @settings[:parallel]
		else
			ri_scf_core =  @settings[:mem] * @settings[:mem_core] / 100 / @settings[:parallel]
		end
		if @settings.set?(:mem_core_correlation)
			ri_corr_core = @settings[:mem] * @settings[:mem_core_correlation] / 100 / @settings[:parallel]
		else
			ri_corr_core = @settings[:mem] * @settings[:mem_core] / 100 / @settings[:parallel]
		end
		# Only for SHM ricc2 calculation, total memory is provided
		if @settings[:parallel_mode] == :shm || @settings[:parallel_mode] == :default
			ri_corr_core *= @settings[:parallel]
		end


		# Other methods than HF: Define is run second time to update the control with method-specific setup
		# For some reason, it does not work when combined into one input
		
		# Load input templates
		method_templates = []
		case @settings[:method]
		when :hf
			if @scf_ri
				# RI-HF setup
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_rihf.erb"))
			else
				# Normal HF does not need additional setup
			end
		when :dft
			if @scf_ri
				# RI-DFT
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_ridft.erb"))
			else
				# Normal DFT
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_dft.erb"))
			end

			# Double-hybrid functionals
			if @settings[:functional] == "b2-plyp"
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_custom_cbas.erb")) if @settings.set?(:auxiliary_basis_mp2)
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_ricc2_rimp2.erb"))
			end

		when :mp2
			# Only RI-MP2 enabled, checked above
			if @scf_ri
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_rihf.erb"))
			end
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_custom_cbas.erb")) if @settings.set?(:auxiliary_basis_mp2)
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_ricc2_rimp2.erb"))
			if @settings[:explicit_correlation] != :none
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_f12.erb"))
			end
		when :mp3
			if @scf_ri
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_rihf.erb"))
			end
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_custom_cbas.erb")) if @settings.set?(:auxiliary_basis_mp2)
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_ricc2_mp3.erb"))
		when :"ccsd"
			if @scf_ri
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_rihf.erb"))
			end
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_custom_cbas.erb")) if @settings.set?(:auxiliary_basis_mp2)
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_ricc2_ccsd.erb"))
		when :"ccsd(t)"
			if @scf_ri
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_rihf.erb"))
			end
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_custom_cbas.erb")) if @settings.set?(:auxiliary_basis_mp2)
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_ricc2_ccsdt.erb"))
		when :rpa
			# DFT setup
			if @scf_ri
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_ridft.erb"))
			else
				method_templates << ERB.new(IO.read(interface_dir + "/templates/define_dft.erb"))
			end
			# RPA setup
			method_templates << ERB.new(IO.read(interface_dir + "/templates/define_rirpa.erb"))
		end

		# DFT-specific setup
		if [:dft, :rpa].include?(@settings[:method])
			require "classes/calculation/dft_setup.rb"
			dft_setup = DftSetup.new(interface_dir + "/dft_functionals.yaml", interface_dir + "/dft_grids.yaml")
		end

		method_templates.each_index{|i|
			# Create input for define
			File.open(calc_dir+"/define2_#{i}.in","w+") {|file|
				file.print(method_templates[i].result(binding))
			}
			# Run define
			system("cd #{calc_dir}; export TURBODIR=#{@settings[:turbomole_turbodir]}; cat define2_#{i}.in | #{@settings[:turbomole_bin_dir]}/define > define2_#{i}.out 2>&1");
			# Look for errors in define output
			unless Grep::file_unicode(calc_dir+"/define2_#{i}.out", /define ended normally/)
				Cuby::error("Turbomole define script failed, see #{calc_dir}/define2_#{i}.out for details")
			end
		}

		# Electrostatic filed
		if @settings.set?(:electric_field)
			template_ef = ERB.new(IO.read(interface_dir + "/templates/define_electric_field.erb"))
			File.open(calc_dir+"/define_ef.in","w+") {|file| file.print(template_ef.result(binding))}
			system("cd #{calc_dir}; export TURBODIR=#{@settings[:turbomole_turbodir]}; cat define_ef.in | #{@settings[:turbomole_bin_dir]}/define > define_ef.out 2>&1");
			Cuby::error("Turbomole define script failed, see #{calc_dir}/define_ef.out for details") unless Grep::file_unicode(calc_dir+"/define_ef.out", /define ended normally/)
		end

		# Properties
		if @settings[:properties].include?(:static_polarizability)
			template_ef = ERB.new(IO.read(interface_dir + "/templates/define_properties_polarizability.erb"))
			fn_part = "prop_polari"
			File.open(calc_dir+"/define_#{fn_part}.in","w+") {|file| file.print(template_ef.result(binding))}
			system("cd #{calc_dir}; export TURBODIR=#{@settings[:turbomole_turbodir]}; cat define_#{fn_part}.in | #{@settings[:turbomole_bin_dir]}/define > define_#{fn_part}.out 2>&1");
			Cuby::error("Turbomole define script failed, see #{calc_dir}/define_#{fn_part}.out for details") unless Grep::file_unicode(calc_dir+"/define_#{fn_part}.out", /define ended normally/)
		end
	end

	def turbomole_cosmoprep
		# Load ERB template
		template = ERB.new(IO.read(interface_dir + "/templates/cosmoprep.erb"))
		# Parse it and write to file
		File.open(calc_dir+"/cosmoprep.in","w+") {|file|
			file.print(template.result(binding))
		}
		# Run cosmoprep
		system("cd #{calc_dir}; export TURBODIR=#{@settings[:turbomole_turbodir]}; cat cosmoprep.in | #{@settings[:turbomole_bin_dir]}/cosmoprep > cosmoprep.out 2>&1");
		# Look for errors in cosmoprep output
		#--------------------------------------------------
		# unless Grep::file_unicode(calc_dir+"/define2_#{i}.out", /define ended normally/)
		# 	Cuby::error("Turbomole define script failed, see #{calc_dir}/define2_#{i}.out for details")
		# end
		#-------------------------------------------------- 
	end

	def turbomole_fix_control
		control = IO.readlines(calc_dir+"/control")
		f = File.open(calc_dir+"/control", "w+")
		control.each{|line|
			if line =~ /^\$ricc2/ && @settings[:spin_component_scaling]
				# SCS caling for RICC2 module
				Cuby::error "Turbomole: spin_component_scaling and sos_scaling can not be set at the same time" if @settings[:sos_scaling]
				f.print line
				f.puts "  scs  cos=#{@settings[:scs_os]}d0  css=#{@settings[:scs_ss]}d0"
			elsif line =~ /^\$ricc2/ && @settings[:sos_scaling]
				# SOS caling for RICC2 module
				Cuby::error "Turbomole: spin_component_scaling and sos_scaling can not be set at the same time" if @settings[:spin_component_scaling]
				# Laplace transformation setup first
				if @settings[:sos_scaling_laplace] != 0
					f.puts "$laplace"
					f.puts "  conv=#{@settings[:sos_scaling_laplace]}"
				end
				# Then the ricc2 keyword
				f.print line
				f.puts "  sos"
			elsif line =~ /^\$scforbitalshift/ && @settings.set?(:scf_levelshift)
				f.puts "$scforbitalshift  automatic=#{@settings[:scf_levelshift]}"
			elsif line =~ /^\$ricc2/ && @settings[:method] == :"ccsd" && @settings[:turbomole_version] <= 6.2
				# Turbomole 6.2 define does not add the "ccsd" line to control
				f.print line
				f.puts "  ccsd"
			elsif line =~ /^\$ricc2/ && @settings[:method] == :"ccsd(t)" && @settings[:turbomole_version] <= 6.2
				# Turbomole 6.2 define does not add the "ccsd(t)" line to control
				f.print line
				f.puts "  ccsd(t)"
			elsif line =~ /^\$drvopt/
				f.print line
				# Switch on calculation of gradient on point charges
				if @point_charges
					f.puts "   point charges"
				end
			elsif line =~ /^\$end/
				# Add things to the end of control
				if @what.include?(:hessian)
					# Hessian calculation - force output to files (automatic in version >= 6.3)
					f.puts "$hessian (projected)   file=hessian"
					f.puts "$vibrational normal modes   file=vib_normal_modes"
					f.puts "$vibrational spectrum   file=vibspectrum"
					f.puts "$dipgrad   file=dipgrad"

					# Hessian calculation - save unprojected hessian
					f.puts "$noproj"
					f.puts "$nprhessian   file=hessian_npr"
					f.puts "$nprvibrational normal modes   file=vib_normal_modes_npr"
					f.puts "$nprvibrational spectrum   file=vibspectrum_npr"
				end
				# Atomic charges
				if what.include?(:atomic_charges)
					case @settings[:atomic_charges]
					when :mulliken
						f.puts "$pop"
					when :nbo
						f.puts "$pop nbo"
					end
				end
				if @point_charges
					f.puts "$point_charges file=point_charges"
					f.puts "$point_charge_gradients file=point_charge_gradients"
				end

				# properties
				if @settings[:properties].include?(:static_polarizability)
					f.puts "$electronic polarizability   file=#{FILE_POLARIZABILITY}"
				end

				# All orbitals instead of just some
				if @what.include?(:mos)
					f.puts "$moprint" 
				end
				
				# Add user-defined keywords
				f.puts @settings[:turbomole_keywords] if @settings.set?(:turbomole_keywords)
				f.print line
			else
				f.print line
			end
		}
		f.close
	end

	def turbomole_run_calculation
		# Write current geometry
		GeoFile::TurbomoleCoord.write(@geometry, :file => calc_dir+"/coord")

		# Write point charges
		if @point_charges
			write_point_charges
		end

		# Get list of executables to be run for selected method
		run_modules = []
		extra_scf_gradient = ''
		if [:hf, :dft, :rpa].include?(@settings[:method])
			if @scf_ri
				run_modules = ["ridft"]
				run_modules << "rdgrad" if @what.include?(:gradient)
			else
				run_modules = ["dscf"]
				run_modules << "grad" if @what.include?(:gradient)
			end

			if @settings[:method] == :dft && @settings[:functional] == "b2-plyp"
				run_modules << "ricc2"
			end

			if @settings[:method] == :rpa
				run_modules << "rirpa"
			end
		elsif [:mp2, :mp3, :ccsd, :"ccsd(t)"].include?(@settings[:method])
			# (RI-MP2 only)
			if @scf_ri
				run_modules = ["ridft"]
				if @what.include?(:gradient) && @settings[:turbomole_scf_grad]
					run_modules << "rdgrad"
					extra_scf_gradient = "rdgrad"
				end

			else
				run_modules = ["dscf"]
				if @what.include?(:gradient) && @settings[:turbomole_scf_grad]
					run_modules << "grad"
					extra_scf_gradient = "grad"
				end
			end
			run_modules << "ricc2"
		end

		if @what.include?(:hessian)
			run_modules << "aoforce"
		end

		# Properties
		run_modules << "escf" if @settings[:properties].include?(:static_polarizability)

		# Delete files with results so that new ones are created
		FileUtils.rm_f(calc_dir+"/energy")
		FileUtils.rm_f(calc_dir+"/gradient")
		FileUtils.rm_f(calc_dir+"/scf_gradient") if @settings[:turbomole_scf_grad]


		# Run calculations
		run_modules.each{|mod|
			# Parallel setup
			if @settings[:parallel] == 1 
				bin_dir = @settings[:turbomole_bin_dir]
				parallel = ""
			elsif @settings[:serial_scf] && ["dscf","ridft","rdgrad","grad"].include?(mod)
				# SCF will be serial / multithreaded only
				bin_dir = @settings[:turbomole_bin_dir]
				parallel = "export PARNODES=#{@settings[:parallel]}; "
			else
				mode = @settings[:parallel_mode]
				mode = :shm if mode == :default
				case mode
				when :shm
					bin_dir = @settings[:turbomole_bin_dir_smp]
					parallel = "export PARNODES=#{@settings[:parallel]}; export PARA_ARCH=SMP;"
				when :mpi
					bin_dir = @settings[:turbomole_bin_dir_mpi]
					parallel = "export PARNODES=#{@settings[:parallel]}; export PARA_ARCH=MPI;"
				when :ga
					bin_dir = @settings[:turbomole_bin_dir_ga]
					parallel = "export PARNODES=#{@settings[:parallel]}; export PARA_ARCH=MPI;"
				end
			end

			system("cd #{calc_dir}; #{parallel} export PATH=#{@settings[:turbomole_scripts_dir]}:$PATH; export TURBODIR=#{@settings[:turbomole_turbodir]}; #{bin_dir}/#{mod} > #{mod}.out 2>&1");
			if mod == extra_scf_gradient
				system("cd #{calc_dir}; mv gradient scf_gradient")
			end
			# Check for errors
			Cuby::error "SCF not converged (calculation #{@name})" if FileTest.exist?(calc_dir+"/dscf_problem")

			# Check the output for errors
		}

	end

	def write_point_charges
		f = File.open(in_calc_dir("point_charges"), "w+")
		f.puts "$point_charges"
		@point_charges.each { |pcharge|
			f.puts "   #{pcharge.x * ANGSTROM2BOHR} #{pcharge.y * ANGSTROM2BOHR} #{pcharge.z * ANGSTROM2BOHR} #{pcharge.charge}\n"	
		}
		f.puts "$end"
		f.close
	end

	def turbomole_read_gradient(filename)
		gradient = Gradient.new
		File.open(filename,"r"){|gradfile|
			# 2 lines header
			2.times {gradfile.gets}
			# Coordinates
			@geometry.size.times {gradfile.gets}
			# Gradient
			@geometry.size.times {
				g = gradfile.gets.strip.gsub("D","e").split.map{|x| x.to_f}
				gradient << Coordinate[g[0],g[1],g[2]] * HARTREE2KCAL / BOHR2ANGSTROM
			}
		}
		return gradient
	end

	def turbomole_read_point_charges_gradient(filename)
		gradient = Gradient.new
		File.open(filename,"r"){|gradfile|
			# 1 line header
			gradfile.gets
			# Gradient
			@point_charges.size.times {
				g = gradfile.gets.strip.gsub("D","e").split.map{|x| x.to_f}
				gradient << Coordinate[g[0],g[1],g[2]] * HARTREE2KCAL / BOHR2ANGSTROM
			}
		}
		return gradient
	end

	def turbomole_read_hessian_from_file(filename, projected = true)
		rows = []
		(@geometry.size * 3).times {|i| rows[i] = []}
		f = File.open(in_calc_dir(filename))
		read = false
		f.each_line{|line|
			if line =~ /\$hessian \(projected\)/ && projected || line =~ /\$nprhessian/ && !projected
				read = true
				next
			end
			if line =~ /^\$/ && read
				break
			end
			if read
				rowindex = line[0..2].to_i - 1
				line_a = line[5..-1].chomp.split
				line_a.each_index {|i|
					rows[rowindex] << line_a[i].to_f * HARTREE2KCAL / BOHR2ANGSTROM / BOHR2ANGSTROM
				}

			end

		}
		f.close
		hessian = Hessian.rows(rows)
		return hessian
	end

	def read_atomic_charges(filename)
		parser = OutputParser.new(filename, @name)
		charges = AtomicCharges.new(@settings[:atomic_charges])
		case @settings[:atomic_charges]
		when :nbo
			parser.add_block(
				:atomic_charges,
				/atomic populations from total density/, 
				/^ *[0-9]+ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
				/=====/,
				{:get => :submatch, :type => :float}
			)
		when :mulliken
			parser.add_block(
				:atomic_charges,
				/atomic populations from total density/, 
				/^ *[0-9][^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
				/moments/,
				{:get => :submatch, :type => :float}
			)
		end
		parser.execute{|name, result, count|
			charges << result
		}
		return charges
	end

	def turbomole_read_polarizability
		lines = IO.readlines(in_calc_dir(FILE_POLARIZABILITY)).map!{|s| s.strip.split}
		polarizability = Polarizability.new("a.u.")

		polarizability.matrix[0,0] = lines[1][0].to_f

		polarizability.matrix[1,0] = polarizability.matrix[0,1] = lines[2][0].to_f
		polarizability.matrix[1,1] = lines[2][1].to_f

		polarizability.matrix[2,0] = polarizability.matrix[0,2] = lines[3][0].to_f
		polarizability.matrix[2,1] = polarizability.matrix[1,2] = lines[3][1].to_f
		polarizability.matrix[2,2] = lines[3][2].to_f

		return polarizability
	end

	def turbomole_read_dipole(results, tm_module)
		case tm_module
		when :ridft
			filename = in_calc_dir("ridft.out")
		when :dscf
			filename = in_calc_dir("dscf.out")
		end

		if filename
			parser = OutputParser.new(filename, @name)
			dipole = []
			parser.add_block(
				:dipole,
				/^ *dipole moment *$/, 
				/^ *[x-z] +[-]*[0-9]*\.[0-9]+ +[-]*[0-9]*\.[0-9]+ +([-]*[0-9]*\.[0-9]+)/, 
				/debye/,
				{:get => :submatch, :type => :float}
			)
			parser.execute {|name, result, count|
				case name
				when :dipole
					dipole[count] = result
				end
			}
			results.multipoles[:dipole] = Multipoles::Dipole.new(
				dipole[0] * AU2DEBYE * DEBYE2CUBY,
				dipole[1] * AU2DEBYE * DEBYE2CUBY,
				dipole[2] * AU2DEBYE * DEBYE2CUBY
			)
		end
		return nil
	end

	def turbomole_read_orbitals(results, tm_module)
		results.molecular_orbitals = MolecularOrbitals.new

		case tm_module
		when :ridft
			filename = in_calc_dir("ridft.out")
		when :dscf
			filename = in_calc_dir("dscf.out")
		end

		if filename
			parser = OutputParser.new(filename, @name)
			moblocks = []
			parser.add_block(
				:mos,
				/^ *orbitals \$scfmo/, 
				/\S/, 
				/=====/,
				{:get => :line}
			)
			this_block = nil
			parser.execute {|name, result, count|
				case name
				when :mos
					if result =~ /irrep/
						# Initialize new block
						this_block = []
						moblocks << this_block
						# Read the line
						values = result.strip.split
						values.shift # remove line name
						this_block << values
					end

					if result =~ /eigenvalues H/
						# Read the line
						values = result.strip.split
						values.shift # remove line name
						values.shift # remove line name
						values.map!{|x| x.to_f}
						this_block << values

					end

					if result =~ /occupation/
						# Read the line
						values = result.strip.split
						values.shift # remove line name
						values.map!{|x| x.to_f}
						this_block << values
					end
				end
			}

			moblocks.each{|block|
				names, energies, occupations = block
				names.each_index{|i|
					name = names[i]
					energy = energies[i]
					if occupations && occupations[i]
						occ = occupations[i]
					else
						occ = 0.0
					end

					results.molecular_orbitals << MolecularOrbital.new(name, energy, occ)
				}
			}
		end
		results.molecular_orbitals.sort_by_e!
		return nil
	end

	def turbomole_read_results
		results = Results.new

		# Read the energy file
		energy_file_data = nil
		IO.readlines(calc_dir+"/energy").each {|line|
			if line =~ / *[0-9]/
				energy_file_data = line.strip.split
			end
		}

		# Read the gradient file and save the gradient
		if @what.include?(:gradient)
			results.gradient = turbomole_read_gradient(in_calc_dir("gradient"))
			if @settings[:turbomole_scf_grad]
				results.gradient_components[:scf_gradient] = turbomole_read_gradient(in_calc_dir("scf_gradient"))
				results.gradient_components[:correlation_gradient] = results.gradient - results.gradient_components[:scf_gradient]
			end

			# Read gradient on point charges
			if @point_charges
				@point_charges.gradient = turbomole_read_point_charges_gradient(in_calc_dir("point_charge_gradients"))
			end
		end

		# Some results must be read from the ricc2 output file
		if @settings[:method] == :dft && @settings[:functional] == "b2-plyp"
			ricc_reader = OutputParser.new(calc_dir+"/ricc2.out", @name)
			ricc_reader.add_pattern(:energy_b2plyp, /^ +\* *B2-PLYP energy.*: +(\S+)/, {:get => :match_array})
			ricc_reader.execute
		end
		if @settings[:method] == :mp3
			ricc_reader = OutputParser.new(calc_dir+"/ricc2.out", @name)
			ricc_reader.add_pattern(:corr_mp2, /^ +\* *MP2 correlation energy.*: +(\S+)/, {:get => :match_array, :multi => :last})
			ricc_reader.add_pattern(:corr_mp2_os_ss, /^.*E\(OS\) += +(\S+) +E\(SS\) += +(\S+)/, {:get => :match_array})
			ricc_reader.add_pattern(:corr_mp3, /^ +\* *MP3 correlation energy.*: +(\S+)/, {:get => :match_array})
			ricc_reader.execute
		end

		if @settings[:method] == :"ccsd"
			ricc_reader = OutputParser.new(calc_dir+"/ricc2.out", @name)
			ricc_reader.add_pattern(:corr_mp2, /^ +\* +MP2 correlation energy.*: +(\S+)/, {:get => :match_array})
			ricc_reader.add_pattern(:corr_mp2_os_ss, /^.*E\(OS\) += +(\S+) +E\(SS\) += +(\S+)/, {:get => :match_array})
			ricc_reader.add_pattern(:corr_ccsd, /^ +\* +correlation energy +: +(\S+)/, {:get => :match_array})
			ricc_reader.execute
		end

		if @settings[:method] == :"ccsd(t)"
			ricc_reader = OutputParser.new(calc_dir+"/ricc2.out", @name)
			ricc_reader.add_pattern(:corr_mp2, /^ +\* +MP2 correlation energy.*: +(\S+)/, {:get => :match_array})
			ricc_reader.add_pattern(:corr_mp2_os_ss, /^.*E\(OS\) += +(\S+) +E\(SS\) += +(\S+)/, {:get => :match_array})
			ricc_reader.add_pattern(:corr_ccsd, /^ +\* +CCSD correlation energy +: +(\S+)/, {:get => :match_array})
			ricc_reader.add_pattern(:corr_ccsdt, /^ +\* +total correlation energy +: +(\S+)/, {:get => :match_array})
			ricc_reader.add_pattern(:ccsdt_e5, /^ +\* +E5 singles and triples contribution +: +(\S+)/, {:get => :match_array})
			ricc_reader.execute
		end

		if @settings[:method] == :mp2
			ricc_reader = OutputParser.new(calc_dir+"/ricc2.out", @name)
			ricc_reader.add_pattern(:corr_mp2, /^ +\* +MP2 correlation energy.*: +(\S+)/, {:get => :match_array, :required => false})
			ricc_reader.add_pattern(:corr_total, /^ +\* +total correlation energy.*: +(\S+)/, {:get => :match_array, :required => false})
			# ricc_reader.add_pattern(:total_e_mp2, /^ +Total Energy +: +(\S+)/, {:get => :match_array, :required => false})
			if @settings[:explicit_correlation] != :none
				ricc_reader.add_pattern(:f12_hf_correction, / +\* +correction to HF from CABS singles.*: +(\S+)/, {:get => :match_array, :required => false})
				ricc_reader.add_pattern(:f12_mp2_correction, / +\* +fixed.*contr\..*: +(\S+)/, {:get => :match_array, :required => false})
			end
			if @settings[:spin_component_scaling] && !(@what.include?(:gradient))
				ricc_reader.add_pattern(:corr_mp2_os_ss, /^.*E\(OS\) += +(\S+) +E\(SS\) += +(\S+)/, {:get => :match_array, :required => false})
				if @settings[:explicit_correlation] != :none
					ricc_reader.add_pattern(:corr_mp2f12_os_ss, /^.*F12\(OS\) += +(\S+) +F12\(SS\) += +(\S+)/, {:get => :match_array, :required => false})
				end
			end
			ricc_reader.execute
		end

		# Save the energy
		case @settings[:method]
		when :hf
			results.energy = energy_file_data[1].to_f * HARTREE2KCAL
			results.atomic_charges = read_atomic_charges(@scf_ri ? calc_dir+"/ridft.out" : calc_dir+"/dscf.out") if @what.include?(:atomic_charges)
		when :dft
			results.energy = energy_file_data[1].to_f * HARTREE2KCAL
			results.atomic_charges = read_atomic_charges(@scf_ri ? calc_dir+"/ridft.out" : calc_dir+"/dscf.out") if @what.include?(:atomic_charges)
			if @settings[:functional] == "b2-plyp"
				results.energy = ricc_reader[:energy_b2plyp][0] * HARTREE2KCAL
			end
		when :rpa
			results.energy = energy_file_data[1].to_f * HARTREE2KCAL
			results.energy_components[:scf_energy] = results.energy
			e_corr_rpa = energy_file_data[4].to_f * HARTREE2KCAL
			results.energy_components[:correlation_rpa] = e_corr_rpa
			results.energy += e_corr_rpa
		when :mp2
			# SCF energy
			results.energy = energy_file_data[1].to_f * HARTREE2KCAL
			results.energy_components[:scf_energy] = results.energy
			# MP2 energy
			begin
				ec_mp2 = ricc_reader[:corr_total][0] * HARTREE2KCAL
			rescue
				begin
					ec_mp2 = ricc_reader[:corr_mp2][0] * HARTREE2KCAL
				rescue
					# Last option: read it from the energy file
					ec_mp2 = energy_file_data[4].to_f * HARTREE2KCAL
				end

			end
			results.energy += ec_mp2

			# Plain MP2 or scaled variants: name the component accordingly
			if @settings[:spin_component_scaling]
				results.energy_components[:correlation_scs_mp2] = ec_mp2
			elsif @settings[:sos_scaling]
				results.energy_components[:correlation_sos_mp2] = ec_mp2
			else
				results.energy_components[:correlation_mp2] = ec_mp2
			end

			# SCS-MP2 - spin components, available only in energy calculation
			# Gradient calcualtion yields only the final result
			if !(@what.include?(:gradient)) && @settings[:spin_component_scaling] && ricc_reader[:corr_mp2_os_ss]
				if ricc_reader[:corr_mp2f12_os_ss]
					results.energy_components[:correlation_mp2_f12_os] = ricc_reader[:corr_mp2f12_os_ss][0] * HARTREE2KCAL
					results.energy_components[:correlation_mp2_f12_ss] = ricc_reader[:corr_mp2f12_os_ss][1] * HARTREE2KCAL
					results.energy_components[:correlation_mp2_f12] = results.energy_components[:correlation_mp2_f12_os] + results.energy_components[:correlation_mp2_f12_ss]
				end
				results.energy_components[:correlation_mp2_os] = ricc_reader[:corr_mp2_os_ss][0] * HARTREE2KCAL
				results.energy_components[:correlation_mp2_ss] = ricc_reader[:corr_mp2_os_ss][1] * HARTREE2KCAL
				results.energy_components[:correlation_mp2] = results.energy_components[:correlation_mp2_os] + results.energy_components[:correlation_mp2_ss]
			end

			# F12
			if @settings[:explicit_correlation] != :none
				results.energy_components[:f12_hf_correction] = ricc_reader[:f12_hf_correction][0] * HARTREE2KCAL
				results.energy_components[:f12_mp2_correction] = ricc_reader[:f12_mp2_correction][0] * HARTREE2KCAL
				results.energy_components[:scf_energy_f12] = results.energy_components[:scf_energy] + results.energy_components[:f12_hf_correction]
				results.energy_components[:correlation_mp2_f12] = results.energy_components[:correlation_mp2] 
				results.energy_components[:correlation_mp2] -= results.energy_components[:f12_mp2_correction]
				results.energy += results.energy_components[:f12_hf_correction]
			end

		when :mp3
			results.energy_components[:scf_energy] = energy_file_data[1].to_f * HARTREE2KCAL
			# MP2
			results.energy_components[:correlation_mp2] = ricc_reader[:corr_mp2][0] * HARTREE2KCAL
			# MP2 - spin components
			results.energy_components[:correlation_mp2_os] = ricc_reader[:corr_mp2_os_ss][0] * HARTREE2KCAL
			results.energy_components[:correlation_mp2_ss] = ricc_reader[:corr_mp2_os_ss][1] * HARTREE2KCAL
			# MP3
			results.energy_components[:correlation_mp3] = ricc_reader[:corr_mp3][0] * HARTREE2KCAL + results.energy_components[:correlation_mp2] 
			# Total
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:correlation_mp3]
		when :"ccsd"
			# SCF energy
			results.energy_components[:scf_energy] = energy_file_data[1].to_f * HARTREE2KCAL
			# MP2 - same as MP3
			results.energy_components[:correlation_mp2] = ricc_reader[:corr_mp2][0] * HARTREE2KCAL
			# MP2 - spin components - same as MP3
			results.energy_components[:correlation_mp2_os] = ricc_reader[:corr_mp2_os_ss][0] * HARTREE2KCAL
			results.energy_components[:correlation_mp2_ss] = ricc_reader[:corr_mp2_os_ss][1] * HARTREE2KCAL
			# CCSD
			results.energy_components[:correlation_ccsd] = ricc_reader[:corr_ccsd][0] * HARTREE2KCAL
			# Total
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:correlation_ccsd]
		when :"ccsd(t)"
			# SCF energy
			results.energy_components[:scf_energy] = energy_file_data[1].to_f * HARTREE2KCAL
			# MP2 - same as MP3
			results.energy_components[:correlation_mp2] = ricc_reader[:corr_mp2][0] * HARTREE2KCAL
			# MP2 - spin components - same as MP3
			results.energy_components[:correlation_mp2_os] = ricc_reader[:corr_mp2_os_ss][0] * HARTREE2KCAL
			results.energy_components[:correlation_mp2_ss] = ricc_reader[:corr_mp2_os_ss][1] * HARTREE2KCAL
			# CCSD
			results.energy_components[:correlation_ccsd] = ricc_reader[:corr_ccsd][0] * HARTREE2KCAL
			# CCSD(T)
			ec_ccsdt = ricc_reader[:corr_ccsdt][0] * HARTREE2KCAL
			results.energy_components[:"correlation_ccsd[t]"] = ec_ccsdt - ricc_reader[:ccsdt_e5][0] * HARTREE2KCAL
			results.energy_components[:"correlation_ccsd(t)"] = ec_ccsdt
			# Total
			results.energy = results.energy_components[:scf_energy] + ec_ccsdt
		end

		# read dipole
		if [:hf, :dft].include?(@settings[:method])
			if @scf_ri
				turbomole_read_dipole(results, :ridft)
			else
				turbomole_read_dipole(results, :dscf)
			end
		end

		# read MOs
		if @what.include?(:mos)
			if @scf_ri
				turbomole_read_orbitals(results, :ridft)
			else
				turbomole_read_orbitals(results, :dscf)
			end
		end


		# Read the hessian
		if @what.include?(:hessian)
			#results.hessian = turbomole_read_hessian_from_file("hessian", true) # Projected hessian
			results.hessian = turbomole_read_hessian_from_file("hessian_npr", false) # Unprojected hessian
		end

		# Gradient on point charges
		if @what.include?(:point_charges_gradient)
			results.point_charges_gradient = @point_charges.gradient
		end

		# properties
		if @settings[:properties].include?(:static_polarizability)
			results.static_polarizability = turbomole_read_polarizability
		end

		return results
	end
end
