################################################################################
#
# DFTB interface
#
# Author: Jan Rezac
# Date created: 2012-08-15
# License: Cuby4 license
# Description: Interface for external calculations in DFTB and DFTB+
# Status: Works
#
################################################################################

#===============================================================================
# Interface to the DFTB+ program
# Free for non-profit use, available after registration
# http://www.dftb-plus.info
#===============================================================================

require "erb" # ERB templating system, used for construction of the input
require "classes/tools/output_parser.rb"
require "fileutils"
require "interfaces/dftb/classes/dftb_slko_files.rb"
require "interfaces/dftb/classes/hubbard_derivs.rb"
require "interfaces/dftb/classes/cpe_parameters.rb"

module InterfaceDftb
	#=======================================================================
	# Interface header
	#=======================================================================
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK, more features to be implemented"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :point_charges, :point_charges_gradient, :atomic_charges]
	MODIFIER = false
	DIRECTORY = "DFTB"
	# Methods provided by the interface:
	METHODS = {
		:"scc-dftb" 		=> [:energy, :gradient, :point_charges, :point_charges_gradient, :atomic_charges],
		:"scc-dftb3" 		=> [:energy, :gradient, :point_charges, :point_charges_gradient, :atomic_charges]
	}
	ATOMIC_CHARGES = [:mulliken]
	#=======================================================================
	

	def prepare_interface
		# Check compatibility
		if @settings[:dftb_use_dftbplus]
			# DFTB+ only warnings
		else
			# Old code warnings
			Cuby::error("Diagonal 3rd order term available only in DFTB+") if @settings[:method] == :"scc-dftb3" && !@settings[:dftb_3rd_order_full]
			Cuby::error("Spin-polarized DFTB available only in DFTB+") if @settings[:multiplicity] != 1
			Cuby::warning("Dispersion might not work properly in DFTB3", :once) if @settings[:dftb_dispersion]
		end

		# Dispersion warning
		unless @settings.set?(:dftb_dispersion)
			Cuby::warning("Default value of 'dftb_dispersion' was changed to 'false' on Feb 2, 2016", :once)
		end

		# Initialize the slko files
		@slko_files = DftbSlkoFiles.new(@settings)
		# Check if there are all parameters needed
		@slko_files.check_geometry(@geometry)

		# Load data tables
		dftb_load_data

		# Load Hubbard derivatives
		if @settings[:method] == :"scc-dftb3"
			datafile = interface_dir + "/data/hubbard_derivs.yaml"
			@hubbard_derivs = HubbardDerivs.from_settings(@settings, datafile)
		end

		# Delta-Pauli
		if @settings[:development][:dftb_dp]
			#--------------------------------------------------
			# @dp_exp = 30.0
			# @dp_alpha = {
			# 	:O => 0.54,
			# 	:C => -0.37,
			# 	:H => -0.42
			# }
			#-------------------------------------------------- 

			@dp_exp = @settings[:x0]
			@dp_alpha = {
				:O => @settings[:x1],
				:C => @settings[:x2],
				:H => @settings[:x3]
			}
			@dp_cutoffs = {
				:H => {
					:H => 3.79965
				},
				:C => {
					:H => 4.46826,
					:C => 5.26856
				},
				:O => {
					:H => 4.279885,
					:C => 4.9, #Arbitrary value, will be calculated soon
					:O => 4.61075
				}
			}
		end

		if @settings[:development][:dftb_dp2]
			@dp_xi = @settings[:x0]
			@dp_zeta = @settings[:x1]

			@dp_alpha = {
				:C => @settings[:x2],
				:H => @settings[:x3]
			}
			@dp_cutoffs = {
				:C => {
					:C => @settings[:x4],
					:H => @settings[:x5]
				},
				:H => {
					:H => @settings[:x6]
				}
			}
		end

		# CPE polarization
		if @settings[:dftb_cpe] != :none
			# Read parameter file
			filename = interface_dir + "/data/parameters_cpe.yaml"
			File.open(filename) {|f|
				@cpe_data = YAML::load(f)
			}
			@cpe_data = @cpe_data[@settings[:dftb_cpe]]
		else
			@cpe_data = nil
		end


		if @settings[:dftb_use_dftbplus]
			# Prepare calculation directory or reuse old
			if calc_dir_mkdir("dftb_in.hsd", "dftb.out") == :old
				calc_using_input
				return 
			end
			# Write input
			calc_writing_input
			dftbplus_write_input(calc_dir+"/dftb_in.hsd")
			# Write geometry
			GeoFile::DftbGen.write(@geometry, :file => calc_dir+"/in.gen")
			# Write point charges if applicable
			dftbplus_write_point_charges if @point_charges
			# Compatibility checks
			Cuby::error("CM3 charges not available in DFTB+") if  @settings[:dftb_cm3_charges]
		else
			# Prepare calculation directory or reuse old
			if calc_dir_mkdir("dftb.in", "dftb.out") == :old
				calc_using_input
				return 
			end
			# Write input
			calc_writing_input
			dftbold_write_input(calc_dir+"/dftb.in")
			dftbold_dispersion_input(calc_dir+"/DISPERSION.INP") if @settings[:dftb_dispersion]
			# Write geometry
			GeoFile::DftbGen.write(@geometry, :file => calc_dir+"/in.gen")
			# Write point charges if applicable
			dftbold_write_point_charges if @point_charges
			# Copy the parameter file for CM3 charges
			if  @settings[:dftb_cm3_charges]
				FileUtils.cp(interface_dir + "/MinFile.current", calc_dir + "/MinFile.current")
			end
		end
	end

	def calculate_interface
		if @settings[:dftb_use_dftbplus]
			dftbplus_run_calculation unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/dftb.out")
			results = dftbplus_read_results
		else
			dftbold_run_calculation unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/dftb.out")
			results = dftbold_read_results
		end

		return results
	end

	def cleanup_interface
		calc_dir_delete if @settings[:job_cleanup]
	end

	#=======================================================================
	# Private methods: common
	#=======================================================================
	

	def dftb_load_data
		if @settings.set?(:dftb_data_file)
			filename = File.expand_path(@settings[:dftb_data_file])
		else
			filename = interface_dir + "/dftb_data.yaml"
		end
		File.open(filename) {|f|
			@dftb_data = YAML::load(f)
		}
	end

	#=======================================================================
	# DFTB+
	#=======================================================================

	def dftbplus_write_input(filename)
		# Elements
		elements = @geometry.elements_in_system

		# Use previous charges
		if @settings[:start_from_previous] &&  FileTest.exist?(calc_dir + '/charges.bin')
			read_previous_charges = true
		else
			read_previous_charges = false
		end

		# Load input template
		template = ERB.new(IO.read(interface_dir + "/dftb_in.hsd.erb"), nil, "<>")

		File.open(filename,"w+") {|file|
			file.print(template.result(binding))
		}
	end

	def dftbplus_run_calculation
		# Write current geometry
		GeoFile::DftbGen.write(@geometry, :file => calc_dir+"/in.gen")

		# Write point charges if applicable
		if @point_charges
			dftbplus_write_point_charges
		end

		# Rewrite input when previous charges exist and should be reused
		if !@prev_charges_input_written && @settings[:existing_calc_dir] != :reuse
			calc_writing_input_in_calculate
			dftbplus_write_input(calc_dir+"/dftb_in.hsd")
			@prev_charges_input_written = true
		end

		# Run calculation
		unless system("cd #{calc_dir}; export OMP_NUM_THREADS=#{@settings[:parallel]}; #{@settings[:dftbplus_exe]} > dftb.out 2>&1")
			Cuby::error("DFTB calculation failed", self)
		end
	end

	def dftbplus_read_results
		# Error check in dftb.out
		parser = OutputParser.new(calc_dir+"/dftb.out")
		parser.add_error(/had been ignored by the parser/, "Some input was ignored by DFTB+, use newer version\n(calculation #{@name}, interface DFTB)")
		parser.execute

		# Results are read from detailed.out
		results = Results.new
		parser = OutputParser.new(calc_dir+"/detailed.out")

		# Line containing energy
		parser.add_pattern(:energy, /Total energy:/)
		# Possible errors
		parser.add_error(/SCC is NOT converged/, "SCC not converged (calculation #{@name}, interface DFTB)")

		# Energy components
		parser.add_pattern(:repulsion, /Repulsive energy: *([-]*[0-9]*\.[0-9]+) *H/, {:required => false, :get => :submatch})


		# Dipole
		parser.add_pattern(:dipole,
			/Dipole.*:\s+([-]*[0-9]*\.[0-9]+)\s+([-]*[0-9]*\.[0-9]+)\s+([-]*[0-9]*\.[0-9]+)\s+Debye/,
			:get => :match_array, :type => :float, :required => false
		)

		# Gradient block
		if @what.include?(:gradient)
			results.gradient = Gradient.new
			parser.add_block(:gradient, /Total Forces/, /[0-9]/, /^ *$/i, {:get => :line_words, :type => :float})
		end

		# Gradient on point charges
		if @point_charges && @what.include?(:point_charges_gradient)
			@point_charges.gradient = Gradient.new
			parser.add_block(:pch_gradient, /Forces on external charges/, /[0-9]/, /^ *$/i, {:get => :line_words, :type => :float})
		end

		# Atomic charges
		if @what.include?(:atomic_charges)
			results.atomic_charges = AtomicCharges.new(:mulliken)
			parser.add_block(
				:atomic_charges, 
				/Net atomic charges/, 
				/^\s+[0-9]+\s+([-]*[0-9]*\.[0-9]+)\s+/,
				/^ *$/i, 
				{:get => :submatch, :type => :float}
			)
		end

		parser.execute {|name, result, count|
			# Gradient block
			case name
			when :gradient
				results.gradient << Coordinate.new(result[0], result[1], result[2]) * HARTREE2KCAL / BOHR2ANGSTROM * -1.0
			when :pch_gradient
				@point_charges.gradient << Coordinate.new(result[0], result[1], result[2]) * HARTREE2KCAL / BOHR2ANGSTROM * -1.0
			when :atomic_charges
				results.atomic_charges << result
			end
		}

		# Energy
		results.energy = parser[:energy] * HARTREE2KCAL if  parser[:energy]

		results.energy_components[:repulsion] = parser[:repulsion] * HARTREE2KCAL

		# Dipole
		if parser[:dipole]
			results.multipoles[:dipole] = Multipoles::Dipole.new(
				parser[:dipole][0] * DEBYE2CUBY,
				parser[:dipole][1] * DEBYE2CUBY,
				parser[:dipole][2] * DEBYE2CUBY
			)
		end

		# Gradient on point charges
		if @what.include?(:point_charges_gradient)
			results.point_charges_gradient = @point_charges.gradient
		end

		return results
	end

	def dftbplus_write_point_charges
		File.open(in_calc_dir("extcharges.xyzc"),"w+") { |f|
			@point_charges.each_index {|i|
				pch = @point_charges[i]
				f.printf("%s", "#{pch.x} #{pch.y} #{pch.z} #{pch.charge}")
				f.puts if i < @point_charges.length - 1
			}
		}
	end
	
	#=======================================================================
	# Original DFTB program
	#=======================================================================

	def dftbold_write_input(filename)
		# Elements
		elements = @geometry.elements_in_system

		# Use previous charges
		if @settings[:start_from_previous] &&  FileTest.exist?(calc_dir + '/CHR.DAT')
			read_previous_charges = true
		else
			read_previous_charges = false
		end

		f = File.open(filename, "w+")

		# First line:
		# mode (3 = SD optimization)
		line = "3"
		# max force (for optimization)
		line += " 1.0e-6"
		# scc mode: 1: no SCC, 2: SCC-DFTB, 2nd order, 3: SCC-DFTB, 3rd order S: special
		#--------------------------------------------------
		# if @settings[:method] == :"scc-dftb3"
		# 	line += " 3"
		# else
		# 	line += " 2"
		# end
		#-------------------------------------------------- 
		line += " S"
		# scc convergence
		line += " #{@settings[:dftb_scc_convergence]}"
		# read charges (T / F)
		if read_previous_charges
			line += " T"
		else
			line += " F"
		end
		# dispersion (T / F)
		if @settings[:dftb_dispersion]
			line += " T"
		else
			line += " F"
		end
		# external field
		if @point_charges
			line += " 'CH'"
		else
			line += " 'NO'"
		end
		f.puts line

		# Extra line: special setup (SCC mode, XH gamma,  fudge)
		#-------------------------------------------------- 
		line = ""
		if @settings[:method] == :"scc-dftb3"
			line += "3"
		else
			line += "2"
		end
		if @settings[:dftb_xh_damping]
			line += " T"
		else
			line += " F"
		end
		line += " F" # "fudge"
		#line += " T" # "Serep"
		line += " F" # l-dependent Hubbard
		line += " 0.0" # Spin polarization
		f.puts line

		# Second line: geometry file
		f.puts "'in.gen'"

		# 3rd line: charge
		f.puts @settings[:charge]

		# 4th line: constraints
		f.puts 0

		# 5th line: output geometry
		f.puts "'out.gen'"

		# 6th line: max. momentum (1,2,3...)
		qn_to_num = {"s" => 1, "p" => 2, "d" => 3}
		f.puts elements.map{|e| qn_to_num[@dftb_data.max_momentum[e]] }.join(" ")

		# 7th line: Hubbard derivatives for the elements
		f.puts elements.map{|e| @hubbard_derivs[e] }.join(" ")

		# 8th line: zeta for XH mod (automatically used in DFTB3)
		f.puts @settings[:dftb_xh_damping_exp]

		# SLKO parameters
		elements.each{|ei|
			elements.each{|ej|
				ee = ei.to_s + @slko_files.separator + ej.to_s
				if @slko_files.downcase
					ee.downcase!
				else
					ee.upcase!
				end
				f.puts "'" + @slko_files.path + "/" + ee + @slko_files.suffix + "'"
			}
		}


		# Last line:
		#       v-- stepsize = 0 => singlepoint
		f.puts "0 0.0 #{@settings[:dftb_e_temp]} 0.0 1000"

		f.close
	end

	def dftbold_dispersion_input(filename)
		# Elements
		elements = @geometry.elements_in_system

		f = File.open(filename, "w+")
		elements.each{|e|
			f.puts "#{@dftb_data.dispersion_para[e].polarisations_old_inp}"
		}
		f.close
	end

	def dftbold_run_calculation
		# Write current geometry
		GeoFile::DftbGen.write(@geometry, :file => calc_dir+"/in.gen")

		# Write point charges if applicable
		if @point_charges
			dftbold_write_point_charges
		end

		# Rewrite input when previous charges exist and should be reused
		if !@prev_charges_input_written && @settings[:existing_calc_dir] != :reuse
			calc_writing_input_in_calculate
			dftbold_write_input(calc_dir+"/dftb.in")
			@prev_charges_input_written = true
		end

		# Run calculation
		unless system("cd #{calc_dir}; export OMP_NUM_THREADS=#{@settings[:parallel]}; #{@settings[:dftb_exe]} < dftb.in > dftb.out 2>&1")
			Cuby::error("DFTB calculation failed", self)
		end
	end

	def dftbold_read_results
		results = Results.new

		parser = OutputParser.new(calc_dir+"/dftb.out")

		# Line containing energy
		parser.add_pattern(:energy, /^    1     1 \/ *[0-9]+ *(\S+)/, {:get => :match_array})

		# Possible errors
		# parser.add_error(/SCC is NOT converged/, "SCC not converged (calculation #{@name}, interface DFTB)")

		parser.execute 

		results.energy = parser[:energy][0] * HARTREE2KCAL

		# Gradient from separate file
		if @what.include?(:gradient)
			results.gradient = Gradient.new
			f = File.open(calc_dir+"/FRC.DAT","r")
			while line =f.gets
				if line =~ /^ *[0-9]/
					n,x,y,z = line.strip.split
					results.gradient << Coordinate.new(x.to_f, y.to_f, z.to_f) * HARTREE2KCAL / BOHR2ANGSTROM * -1.0
				end
			end
			f.close
		end

		return results
	end

	def dftbold_write_point_charges
		File.open(in_calc_dir("EXTCHARGES.INP"),"w+") { |f|
			f.puts @point_charges.size
			@point_charges.each_index {|i|
				pch = @point_charges[i]
				f.printf("%s", "#{pch.x} #{pch.y} #{pch.z} #{pch.charge}")
				f.puts if i < @point_charges.length - 1
			}
		}
	end

	#=======================================================================
	# Module classes used to store some data tables
	#=======================================================================
	
	class SpinConstants
		attr_accessor :wss, :wsp, :wsd, :wpp, :wpd, :wdd

		def to_s
			if @wdd
				return "#{@wss} #{@wsp} #{@wsd} #{@wsp} #{@wpp} #{@wpd} #{@wsd} #{@wpd} #{@wdd}"
			elsif @wpp
				return "#{@wss} #{@wsp} #{@wsp} #{@wpp}"
			elsif @wss
				return "#{@wss}"
			else
				raise "Spin constants not defined"
			end
		end
	end

	class DispersionData
		attr_accessor :radius
		attr_accessor :polarisations
		attr_accessor :polarisations_old_inp
	end

	class DFTBData
		attr_accessor :max_momentum # Hash element => string
		attr_accessor :dispersion_para
		attr_accessor :spin_constants
	end

end
