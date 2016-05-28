################################################################################
#
# deMon interface
#
# Author: Jan Rezac
# Date created: 2013-06-16
# License: Cuby4 license
# Description: Interface for external calculations using deMon
# Status: Works
#
################################################################################

#===============================================================================
# Interface to deMon
# http://www.demon-software.com
#===============================================================================

require "classes/calculation/dft_setup.rb"
require "classes/tools/output_parser.rb"
require "fileutils"

module InterfaceDemon
	#=======================================================================
	# Interface header
	#=======================================================================
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "Seems to work"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :atomic_charges, :gradient, :point_charges, :point_charges_gradient, :ghost_atoms]
	MODIFIER = false
	DIRECTORY = "DEMON"
	# Methods provided by the interface:
	METHODS = {
		:"dft"		=> [:energy, :atomic_charges, :gradient, :point_charges, :point_charges_gradient]
	}
	# Atomic charges
	ATOMIC_CHARGES = [:spatial, :mulliken, :loewdin, :hirshfeld_fragment, :becke, :hirshfeld]
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:fragment_N, :optional, "N = 1...demon_fragment_count, blocks defining fragments (charge, multiplicity, selction)"]
	]
	#=======================================================================
	
	# Large files to delete if requested
	DELETE_LARGE_FILES = "deMon.rst deMon.mem iofrg.scr BASIS AUXIS ECPS MCPS FFDS"

	def prepare_interface
		if @what.include?(:gradient) && !@settings[:demon_with_cuby_interface] && @settings[:demon_version] < "4.2.4"
			Cuby::error "Gradient calculation is not supported with the current version of deMon\nwithout a patch adding a cuby interface"
		end

		if @what.include?(:point_charges_gradient) && @settings[:demon_version] < "4.2.4"
			Cuby::error "Calculation of gradient on point charges is available only in demon 4.2.4 and above"
		end

		# Number atoms in geometry
		@atom_numbering = @geometry.atom_numbering_by_element

		# Prepare calculation directory
		if calc_dir_mkdir("deMon.header", "deMon.out") == :old
			calc_using_input # Old input should be reused
			return
		end

		# If basis set is not specified and there is a file corresponding
		# to the basis set in cuby's library, use is (otherwise, the default
		# of the keyword is used later)
		unless @settings.set?(:demon_basis_file)
			# get a list of basis sets in cuby's library
			@cuby_basis_sets = {}
			File.open(interface_dir + "/basissets/basissets.yaml"){|f|
				YAML.load(f).each_pair{|key, value|
					@cuby_basis_sets[key.downcase] = interface_dir + "/basissets/" + value
				}
			}
			# set the keyword if basis set was found
			if @cuby_basis_sets.has_key?(@settings[:basisset].downcase)
				Cuby::log.puts_debug "Demon basis set '#{@settings[:basisset]}' taken from cuby's library"
				@settings[:demon_basis_file] = @cuby_basis_sets[@settings[:basisset].downcase]
			end
		end

		# Copy basis set files
		demon_copy_basis_sets

		# Write the keyword part of the input
		demon_write_keywords # This writes file deMon.header

		# Create complete input by combining deMon.header and the current geometry
		demon_write_input_geo

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		demon_run unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/deMon.out")
		return demon_read
	end

	def cleanup_interface
		if @settings[:job_cleanup]
			# Delete everything
			calc_dir_delete 
		elsif @settings[:delete_large_files]
			# Delete only large files (and some unnecessary files as well)
			system("cd #{calc_dir}; rm -f #{DELETE_LARGE_FILES}")
		end
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def demon_copy_basis_sets
		# Check if the keywords are set
		Cuby::error("Path to the deMon basis file not set (keyword 'demon_basis_file')") if @settings[:demon_basis_file] == ""
		Cuby::error("Path to the deMon auxiliary basis file not set (keyword 'demon_auxis_file')") if @settings[:demon_auxis_file] == ""
                Cuby::error("Path to the deMon effective potential file not set (keyword 'demon_ecps_file')") if @settings[:demon_ecps_file] == ""
                Cuby::error("Path to the deMon model core potential  file not set (keyword 'demon_mcps_file')") if @settings[:demon_mcps_file] == ""

		# All but basis
		# Check if the files exist
		unless FileTest.exist?(@settings[:demon_auxis_file])
			Cuby::error("deMon auxiliary basis set file not found (demon_auxis_file = '#{@settings[:demon_auxis_file]}')")
		end

		# Copy the files to the calculation directory
		FileUtils.cp(@settings[:demon_auxis_file], in_calc_dir("AUXIS"))
                FileUtils.cp(@settings[:demon_ecps_file], in_calc_dir("ECPS"))
                FileUtils.cp(@settings[:demon_mcps_file], in_calc_dir("MCPS"))
		# Create the FFDS file for new demon
		if @settings[:demon_version] >= "4.2.5"
			File.open(in_calc_dir("FFDS"),"w+"){|f| f.puts "Empty file." }
		end

		# Basis set
		if @settings[:demon_basis_file] =~ /;/
			# Multiple files
			Cuby::log.puts_debug("Found multiple basis set files, concatenating them")
			system("rm -rf #{in_calc_dir("BASIS")}")
			files = @settings[:demon_basis_file].strip.split(/\s*;\s*/)
			files.each{|fn|
				unless FileTest.exist?(fn)
					Cuby::error("deMon basis set file not found (demon_basis_file = '#{fn}')")
				end
				system("cat \"#{fn}\" >> #{in_calc_dir("BASIS")}")
			}
		else
			# Check if the files exist
			unless FileTest.exist?(@settings[:demon_basis_file])
				Cuby::error("deMon basis set file not found (demon_basis_file = '#{@settings[:demon_basis_file]}')") 
			end
			# Copy the files to the calculation directory
			FileUtils.cp(@settings[:demon_basis_file], in_calc_dir("BASIS"))
		end
	end

	def demon_write_keywords
		# Elements in system
		elements = @geometry.elements_in_system

		# flags
		use_fragments = false
		if @settings[:demon_fragment_guess]
			use_fragments = true
		end

		# DFT setup
		dft_setup = DftSetup.new(interface_dir + "/dft_functionals.yaml", interface_dir + "/dft_grids.yaml")

		# Write header of the input file
		filename = in_calc_dir("deMon.header")
		f = File.open(filename,"w+")
		
		# Write Title
		f.puts "TITLE demon calculation"

		# Electron configuration
		f.puts "CHARGE #{@settings[:charge]}"
		f.puts "MULTIPLICITY #{@settings[:multiplicity]}"

		# Basis
		f.puts "BASIS (#{@settings[:basisset]})"
		# + elements
		@settings.elements_hash(:basisset_elements).each_pair{|element, basis|
			f.puts "#{element.to_s} (#{basis})" if elements.include?(element)
		}

		# Auxiliary basis
		if @settings.set?(:auxiliary_basis_scf)
			if @settings[:demon_version] >= "4.2.5" &&  @settings.set?(:demon_maxl_auxis)
				  f.puts "AUXIS (#{@settings[:auxiliary_basis_scf]}) LMAX=#{@settings[:demon_maxl_auxis]}"
			else
				  f.puts "AUXIS (#{@settings[:auxiliary_basis_scf]})"
			end
		end
		# + elements
		@settings.elements_hash(:auxiliary_basis_scf_elements).each_pair{|element, basis|
			f.puts "#{element.to_s} (#{basis})" if elements.include?(element)
		}

		# Pseudopotentials
		if @settings.set?(:pseudopotentials)
			f.puts "ECPS (#{@settings[:pseudopotentials]})"
		else
			f.puts "ECPS (NONE)"
		end
		@settings.elements_hash(:pseudopotentials_elements).each_pair{|element, pp|
			f.puts "#{element.to_s} (#{pp})" if elements.include?(element)
		}

		# Integrals
		if @settings[:scf_integral_storage] == :direct
			if  @settings[:scf_multipole_expansion]
				f.puts "ERIS MULTIPOLE DIRECT"
			else
				f.puts "ERIS DIRECT"
			end
		elsif @settings[:scf_integral_storage] == :mixed
			f.puts "ERIS MIXED"
		else # :disc or :default
			f.puts "ERIS CONVENTIONAL"
		end

		# orbitals
		case @settings[:demon_orbitals]
		when :spherical
			f.puts "ORBITAL SPHERICAL"
		when :cartesian
			f.puts "ORBITAL CARTESIAN"
		end
		f.puts "SHIFT #{@settings[:scf_levelshift]}" if @settings.set?(:scf_levelshift)

		# Functional + density fitting
		f.print "VXCTYPE #{dft_setup.functional(@settings)}"
		if [:scf, :both].include?(@settings[:density_fitting])
			f.puts 
		else
			f.puts " BASIS"
		end

		# DFT grid
		if @settings[:dft_grid] == :custom
			f.puts "GRID #{@settings[:dft_grid_custom]}"
		else
			f.puts "GRID #{dft_setup.grid(@settings)}"
		end

		# SCF setup
		if @settings[:spin_restricted]
			scftype = "RKS"
			if @settings[:multiplicity] != 1
				scftype = "ROKS"
			end
		else
			scftype = "UKS"
		end
		f.puts "SCFTYPE #{scftype} TOL=1.0e-#{@settings[:scf_convergence]} CDF=1.0e-#{@settings[:density_convergence]} MAX=#{@settings[:scf_cycles]}"

		# Population analysis
		if @what.include?(:atomic_charges)
			case @settings[:atomic_charges]
			when :spatial
				case @settings[:demon_spatial_population]
				when :plane_divider
					atom1 = GeometrySelections.atomlist_to_string(@geometry.atomlist_from_selection(@settings[:demon_plane_axis_a])).gsub("-",":").gsub(",",";")
					atom2 = GeometrySelections.atomlist_to_string(@geometry.atomlist_from_selection(@settings[:demon_plane_axis_b])).gsub("-",":").gsub(",",";")
					ratio = @settings[:demon_plane_ratio]
					f.puts "POPULA ALSEN=#{atom1}|#{atom2}|#{ratio}"
				end
			when :mulliken
					f.puts "POPULA MULLIKEN"
			when :loewdin
					f.puts "POPULA LOEWDIN"
			when :becke
					f.puts "POPULA BECKE"
			when :hirshfeld
					pop_input = ""
					if  @settings[:demon_version] == "2.5.4"
						  pop_input << "PROA"
					end
					f.puts "POPULA HIRSH #{pop_input}"
			when :hirshfeld_fragment
				f.puts "POPULA HIRSH PROF"
				unless @settings[:demon_fragment_guess]
					Cuby::error "Guess from fragments must be enabled when fragment density is used"
				end
				use_fragments = true
			end
		end

		# Constrained DFT
		cdft_input = ""
		if @settings[:cdft_hda]
			cdft_input << "HDA "
                        cdft_input << "PHASE=#{@settings[:demon_cdft_phase]} " if @settings.set?(:demon_cdft_phase)
		end
                if  @settings[:demon_version] >= "2.5.4"
                        cdft_input << "START=1.0E#{@settings[:demon_cdft_start]} "
                end
		if @settings[:demon_cdft_deformation_density] && @settings[:demon_version] > "2.5.4"
			 cdft_input << "DEFORM "
		end
                minus = ""
                minus = "-" if @settings[:demon_cdft_force_convergence]
		case @settings[:demon_constrained_dft]
		when :spatial
			atom1 = GeometrySelections.atomlist_to_string(@geometry.atomlist_from_selection(@settings[:demon_plane_axis_a])).gsub("-",":").gsub(",",";")
			atom2 = GeometrySelections.atomlist_to_string(@geometry.atomlist_from_selection(@settings[:demon_plane_axis_b])).gsub("-",":").gsub(",",";")
			ratio = @settings[:demon_plane_ratio]

			f.puts "CNSDFT TOL=1.0e-#{@settings[:demon_cdft_tolerance]} ALSEN #{cdft_input}"
			f.puts "0.0 CHARGE ALSEN #{atom1}|#{atom2}|#{ratio}"
		when :mulliken
			f.puts "CNSDFT TOL=#{minus}1.0e-#{@settings[:demon_cdft_tolerance]} MULLIK #{cdft_input}"
			cdft_constraints_input(f)
		when :becke
			f.puts "CNSDFT TOL=#{minus}1.0e-#{@settings[:demon_cdft_tolerance]} BECKE #{cdft_input}"
			cdft_constraints_input(f)
		when :voronoi
			if @settings[:demon_version] <= "2.5.4"
				f.puts "CNSDFT TOL=#{minus}1.0e-#{@settings[:demon_cdft_tolerance]} VDD PROA #{cdft_input}"
				cdft_constraints_input(f)
			else
				f.puts "CNSDFT TOL=#{minus}1.0e-#{@settings[:demon_cdft_tolerance]} VORON #{cdft_input}"
				cdft_constraints_input(f)
			end
		when :loewdin
			f.puts "CNSDFT TOL=#{minus}1.0e-#{@settings[:demon_cdft_tolerance]} LOEWD #{cdft_input}"
			cdft_constraints_input(f)
		when :hirshfeld_fragment
			unless @settings[:demon_fragment_guess]
				Cuby::error "Guess from fragments must be enabled when fragment density is used in CDFT"
			end
			f.puts "CNSDFT TOL=#{minus}1.0e-#{@settings[:demon_cdft_tolerance]} HIRSH PROF #{cdft_input}"
			use_fragments = true
			cdft_constraints_input(f)
		when :hirshfeld_atomic
                        cdft_input << "PROA" if  @settings[:demon_version] == "2.5.4"
			f.puts "CNSDFT TOL=#{minus}1.0e-#{@settings[:demon_cdft_tolerance]} HIRSH #{cdft_input}" 
			cdft_constraints_input(f)
		end

		# Extra keywords
		if @settings.set?(:demon_keywords)
			f.puts @settings[:demon_keywords]
		end

		# Define the fragments
		if use_fragments
			f.puts "FRAGMENT"
			@settings[:demon_fragment_count].times{|i|
				block = @settings.block("fragment_#{i+1}".to_sym)
				f.puts("#{block[:charge]} #{block[:multiplicity]}")
			}
		end

		if @point_charges
			f.puts "EMBED FILE"
			demon_write_pch_file
		end

		# Gradient
		if @what.include?(:gradient)
			# There are multiple options how to get the gradient, starting with the best one:
			if  @settings[:demon_version] >= "4.3.3"
				# Use the QMMM output of the new version - syntax was changed, no CUBY option
				f.puts "QM/MM CHARMM"
			elsif  @settings[:demon_version] >= "4.2.4"
				# Use the QMMM output of the new version
				f.puts "QM/MM#{@settings[:demon_with_cuby_interface] ? ' CUBY' : ''} #{@settings[:demon_cuby_options]}"
			elsif @settings[:demon_with_cuby_interface]
				# Use custom cuby interface patch
				f.puts "CUBY GRAD #{@settings[:demon_cuby_options]}"
			else
				# Get it from optimization - dos not work now, disabled
				f.puts "OPTIMIZATION CARTESIAN MAX=1\nPRINT OPT"
			end
		else
			if @settings[:demon_with_cuby_interface] && @settings[:demon_version] < "4.2.5"
				f.puts "CUBY #{@settings[:demon_cuby_options]}" 
			end
		end



		f.close
	end

	def demon_write_pch_file
		#: Write external point charges in deMon format
		write_element = @settings[:demon_charges_pol]
		if write_element && !@point_charges[0].from_atom
			Cuby::error("Calculation with polarizable external charges requested\nbut polarizability information not found")
		end
		File.open(in_calc_dir("/deMon.cub"),"w+") { |f|
			@point_charges.each_index {|i|
				pch = @point_charges[i]
                                f.printf("%12.5f %12.5f %12.5f %10.5f",pch.x,pch.y,pch.z,pch.charge)
				f.printf("   pol=%s",pch.from_atom.element) if write_element
				f.puts if i < @point_charges.length - 1
			}
		}
		return nil
	end

	CONSTRAINT_ALL_KEYS = ["type", "selection", "target", "donor", "acceptor"]
	CONSTRAINT_REQUIRED_KEYS = ["type", "selection", "target"]
	CONSTRAINT_TARGETS = ["charge", "spin"]

	def cdft_constraints_input(f)
		@settings[:cdft_constraints].each{|cons|
			# Check if there are no invalid keys
			cons.each_key{|k|
				Cuby::error("Invalid constraint property '#{k}'") unless CONSTRAINT_ALL_KEYS.include?(k)
			}
			# Check if all the required keys are there
			CONSTRAINT_REQUIRED_KEYS.each{|k|
				Cuby::error("Incomplete specification of DFT constraint: missing '#{k}'") unless cons.has_key?(k)
			}
			# Check the type
			Cuby::error("Incorrect type of constraint '#{cons['type']}'") unless CONSTRAINT_TARGETS.include?(cons['type'])
			# Target must be a number
			Cuby::error("CDFT: Target must be a number") unless cons['target'].kind_of?(Numeric)

			str = "#{cons['target'].to_f} #{cons['type'].upcase}"
			str << " DO" if cons['donor']
			str << " AC" if cons['acceptor']
			# New selection format
			str << " "
			str << GeometrySelections.atomlist_to_string(@geometry.atomlist_from_selection(cons['selection'].to_s)).gsub("-",":").gsub(",",";")
			f.puts str
		}
	end

	def demon_write_input_geo
		# Write input file using the existing header, only geometry is added
		filename = in_calc_dir("deMon.inp")
		f = File.open(filename,"w+")

		# Header from file
		header = IO.readlines(in_calc_dir("deMon.header"))
		f.puts header

		# Initial guess
		guess_only = ""
		guess_only = " ONLY" if @settings[:demon_no_scf]
		if @settings[:demon_fragment_guess]
			f.puts "GUESS FRAGMENT#{guess_only}"
			use_fragments = true
		else
			if FileTest.exist?(in_calc_dir("deMon.rst")) && @settings[:start_from_previous]
				Cuby.log.puts_debug "Demon starts from a .rst file"
				f.puts "GUESS RESTART#{guess_only}"
			else
				Cuby.log.puts_debug "Demon starts from a new guess"
			end
		end

		# Geometry start
		f.puts "GEOMETRY CARTESIAN ANGSTROM"

		# Fragments
		if @settings[:demon_fragment_guess]
			use_fragments = true
			fragment_info = Array.new(@geometry.size, "")
			#!# Check whether the fragments overlap or not
			@settings[:demon_fragment_count].times{|i|
				selection = @settings.block("fragment_#{i+1}".to_sym)[:selection]
				@geometry.atomlist_from_selection(selection).each{|at_i|
					fragment_info[at_i] = "   F#{i+1}"
				}
			}
		end

		# Geometry
		@geometry.each_index{|i|
			atom = @geometry[i]
		 	f.printf("%2s%-3d %15.10f %15.10f %15.10f", atom.element.to_s, @atom_numbering[i], atom.x, atom.y, atom.z)
			f.printf(fragment_info[i]) if use_fragments
			f.print "   0" if atom.properties[:ghost]
			f.puts
		}

		f.close
	end

	def demon_run
		# Build input with the current geometry
		demon_write_input_geo

		# Write point charges
		demon_write_pch_file if @point_charges

		# Run demon
		command = "cd #{calc_dir};"
		command << "export LD_LIBRARY_PATH=#{@settings[:demon_lib_dir]}:$LD_LIBRARY_PATH;"
		command << "export TMPDIR=`pwd`;"
		if @settings[:parallel] > 1
			# Parallel setup
			if @settings[:demon_call_mpi]
				command << "mpirun -np #{@settings[:parallel]} #{@settings[:demon_exe_mpi]} 2> deMon.err > deMon.stdout;"
			else
				command << "#{@settings[:demon_exe_mpi]} #{@settings[:parallel]} 2> deMon.err > deMon.stdout;"
			end
		else
			# Serial command
			command << "#{@settings[:demon_exe]} deMon.inp 2> deMon.err > deMon.stdout;"
		end
		unless system(command)
			Cuby::error "deMon returned nonzero exit code (calculation #{@name})"
		end
	end

	def demon_read_qmmm_file_gradient
		# Set up parser
		parser = OutputParser.new(in_calc_dir("deMon.qmm"), @name)

		# Gradient
		gradient = Gradient.new
		parser.add_block(
			:gradient,
			/QMFORCES/, 
			/^ +([^ ]*) +([^ ]*) +([^ ]*)/, 
			/^\s*$/,
			{:get => :match_array, :type => :float}
		)

		# Gradient on point charges
		if @point_charges && @what.include?(:point_charges_gradient)
			@point_charges.gradient = Gradient.new
			parser.add_block(
				:pch_gradient,
				/EMBEDFORCES/, 
				/^ +([^ ]*) +([^ ]*) +([^ ]*)/, 
				/^\s*$/,
				{:get => :match_array, :type => :float}
			)
		end

		# Run parser
		parser.execute{|name, result, count|
			case name
			when :gradient
				x,y,z = result
				gradient << Coordinate[x,y,z] * HARTREE2KCAL / BOHR2ANGSTROM
			when :pch_gradient
				x,y,z = result
				@point_charges.gradient << Coordinate[x,y,z] * HARTREE2KCAL / BOHR2ANGSTROM
			end
		}


		return gradient
	end

	def demon_read # => Results

		# Whan a modified demon is used, read the yaml data from output first
		if @settings[:demon_with_cuby_interface]
			yaml_data = demon_read_cuby_data
		end

		# Initialize results
		results = Results.new
		if what.include?(:atomic_charges)
			case @settings[:atomic_charges]
			when :spatial
				results.atomic_charges = AtomicCharges.new(:spatial, @geometry.size)
			when :mulliken
				results.atomic_charges = AtomicCharges.new(:mulliken)
			when :loewdin
				results.atomic_charges = AtomicCharges.new(:loewdin)
			when :becke
				results.atomic_charges = AtomicCharges.new(:becke)
			when :hirshfeld_fragment
				results.atomic_charges = AtomicCharges.new(:hirshfeld)
			when :hirshfeld
				results.atomic_charges = AtomicCharges.new(:hirshfeld)
			end
		end

		# Crate parser
		parser = OutputParser.new(in_calc_dir("deMon.out"), @name)

		# Errors
		parser.add_error(/SCF NOT CONVERGED/, "SCF did not converge")
		parser.add_error(/OPTIMIZATION OF LAGRANGE MULTIPLIER/, "CDFT constraint optimization failed")
		parser.add_error(/NOT IN BASIS SET FILE/, "Basis set not found in the basis set file")
		parser.add_error(/ABNORMAL PROGRAM TERMINATION/, "DeMon calculation eneded with an error")

		# SCF energy
		parser.add_pattern(:e_scf, /^ *TOTAL ENERGY += +([^ ]*)/, :get => :submatch, :multi => :last) unless @settings[:demon_no_scf]

		# Gradient read from output by this parser (only one of the options)
		if @what.include?(:gradient) && @settings[:demon_version] < "4.2.4"
			results.gradient = Gradient.new
			unless @settings[:demon_with_cuby_interface]
				parser.add_block(
					:gradient,
					/CARTESIAN GRADIENT/, 
					/^ +[^ ]+ +[^ ]+ +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+)/, 
					/RMSQ/,
					{:get => :match_array, :type => :float}
				)
			end
		end
		
		# Parser for atomic charges
		if what.include?(:atomic_charges)
			case @settings[:atomic_charges]
			when :spatial
				case @settings[:demon_spatial_population]
				when :plane_divider
					parser.add_pattern(:space_pop_l, /^CHARGE OF THE LEFT REGION *: +([^ ]*)/, :get => :submatch)
					parser.add_pattern(:space_pop_r, /^CHARGE OF THE RIGHT REGION *: +([^ ]*)/, :get => :submatch)
				end
			when :mulliken
				if @settings[:demon_version] =~ /^2\./
					parser.add_block(
						:atomic_charges,
						/MULLIKEN POPULATION ANALYSIS/, 
						/^ +[0-9]+ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
						/ATOMIC SPIN CHARGES/,
						{:get => :submatch, :type => :float}
					)
				else
					parser.add_block(
						:atomic_charges,
						/MULLIKEN POPULATION ANALYSIS/, 
						/^ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
						/ATOMIC VALENCE MATRIX/,
						{:get => :submatch, :type => :float}
					)
				end
			when :loewdin
				parser.add_block(
					:atomic_charges,
					/LOEWDIN POPULATION ANALYSIS/, 
					/^ +[0-9]+ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
					/ATOMIC SPIN CHARGES/,
					{:get => :submatch, :type => :float}
				)
			when :becke
				if @settings[:demon_version] =~ /^2\./
					parser.add_block(
						:atomic_charges,
						/BECKE POPULATION ANALYSIS/, 
						/^ +[0-9]+ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
						/ATOMIC SPIN CHARGES/,
						{:get => :submatch, :type => :float}
					)
				else
					parser.add_block(
						:atomic_charges,
						/BECKE POPULATION ANALYSIS/, 
						/^ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
						/NUMERICAL/,
						{:get => :submatch, :type => :float}
					)
				end
			when :hirshfeld_fragment
				parser.add_block(
					:atomic_charges,
					/HIRSHFELD POPULATION ANALYSIS/, 
					/^ +[0-9]+ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
					/ATOMIC SPIN CHARGES/,
					{:get => :submatch, :type => :float}
				)
			when :hirshfeld
				if @settings[:demon_version] =~ /^2\./
					parser.add_block(
						:atomic_charges,
						/HIRSHFELD POPULATION ANALYSIS/, 
						/^ +[0-9]+ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
						/ATOMIC SPIN CHARGES/,
						{:get => :submatch, :type => :float}
					)
				else
					parser.add_block(
						:atomic_charges,
						/HIRSHFELD POPULATION ANALYSIS/, 
						/^ +[^ ]+ +([-]*[0-9]*\.[0-9]+)/, 
						/NUMERICAL/,
						{:get => :submatch, :type => :float}
					)
				end
			end
		end

		# Run parser
		parser.execute{|name, result, count|
			case name
			when :gradient
				x,y,z = result
				results.gradient << Coordinate[x,y,z] * HARTREE2KCAL / BOHR2ANGSTROM
			when :atomic_charges
				results.atomic_charges << result
			end
		}

		if @what.include?(:gradient)
			if  @settings[:demon_version] >= "4.2.4"
				# Use the QMMM output of the new version
				results.gradient = demon_read_qmmm_file_gradient
			elsif @settings[:demon_with_cuby_interface]
				# Use custom cuby interface patch
				yaml_data["gradient_first_state"].each_index{|i|
					results.gradient[i] = Coordinate.from_array(yaml_data['gradient_first_state'][i]) * (HARTREE2KCAL / BOHR2ANGSTROM)
				}
			end
		end

		# SCF energy
		if @settings[:demon_no_scf]
			# SCF not performed, set energy to 0
			results.energy = 0.0
		else
			results.energy = parser[:e_scf] * HARTREE2KCAL 
		end

		# Atomic charges
		if what.include?(:atomic_charges)
			case @settings[:atomic_charges]
			when :spatial
				case @settings[:demon_spatial_population]
				when :plane_divider
					# The charge in the half-space is put on the first atom in the selection
					atom1 = @geometry.atomlist_from_selection(@settings[:demon_plane_axis_a])[0]
					atom2 = @geometry.atomlist_from_selection(@settings[:demon_plane_axis_b])[0]
					# This has to be fixed:
					results.atomic_charges[atom1] = parser[:space_pop_l]
					results.atomic_charges[atom2] = parser[:space_pop_r]
				end
			end
		end

		# Gradient on point charges
		if @what.include?(:point_charges_gradient)
			results.point_charges_gradient = @point_charges.gradient
		end

		# Print additional data from deMon
		if @settings[:demon_with_cuby_interface]
			# Print additional data from deMon
			if yaml_data && yaml_data['print']
				puts "deMon output:"
				yaml_data['print'].each{|row|
					puts "  #{row['name']}: #{row['value']} #{row['unit']}"
				}
				yaml_data['energy_components'].each_pair{|name, value|
					results.energy_components[name.to_sym] = value
				}
			end
		end

		return results
	end

	def demon_read_cuby_data
		s = ""
		read = false
		f = File.new(in_calc_dir("deMon.out"),"r")
		f.each_line {|line|
			read = false if line =~ /CUBY DATA END/
			s += line if read
			read = true if line =~ /CUBY DATA START/
		}
		return YAML.load(s)
	end
end
