################################################################################
#
# PSI4 interface
#
# Author: Jan Rezac
# Date created: 2013-08-22
# License: Cuby4 license
# Description: Interface for external calculations using PSI4
# Status: Works
#
################################################################################

#===============================================================================
# Interface to PSI4
# http://www.psicode.org
#===============================================================================

require "classes/tools/output_parser.rb"

module InterfacePsi4
	#=======================================================================
	# Interface header
	#=======================================================================
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :ghost_atoms]
	MODIFIER = false
	DIRECTORY = "PSI4"
	# Methods provided by the interface:
	METHODS = {
		:"hf"		=> [:energy, :gradient],
		:"dft"		=> [:energy, :gradient],
		:"mp2"		=> [:energy, :gradient],
		:"mp3"		=> [:energy, :gradient],
		:"ccsd"		=> [:energy],
		:"ccsd(t)"	=> [:energy],
		:"sapt"		=> [:energy]
		#--------------------------------------------------
		# :"qcisd"	=> [:energy],
		# :"qcisd(t)"	=> [:energy]
		#-------------------------------------------------- 
	}
	#SOLVENTS = [:cosmo]
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:molecule_a, :optional, "Definition of the first monomer for SAPT calculation (selection, charge, multiplicity)"],
		InputBlock[:molecule_b, :optional, "Definition of the second monomer for SAPT calculation (selection, charge, multiplicity)"]
	]
	#=======================================================================

	def prepare_interface
		# Prepare calculation directory
		if calc_dir_mkdir("psi4.header", "psi4.out") == :old
			calc_using_input # Old input should be reused
			return
		end

		# Write the keyword part of the input
		psi4_write_keywords # This writes file psi4.header

		# Create complete input by combining psi4.header and the current geometry
		psi4_write_input_geo

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		psi4_run unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/psi4.out")
		return psi4_read
	end

	def cleanup_interface
		calc_dir_delete if @settings[:job_cleanup]
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def psi4_write_keywords
		# Write header of the input file (geometry is inserted later)
		filename = in_calc_dir("psi4.header")
		f = File.open(filename,"w+")

		# Memory
		f.puts "memory #{@settings[:mem]} mb"

		# Molecule information
		f.puts "molecule {"
		f.puts "   symmetry c1" unless @settings[:use_symmetry]
		f.puts "   units angstrom"
		f.puts "   %geometry%"
		f.puts "}"

		# RHF/UHF reference
		if @settings[:spin_restricted]
			f.puts "set reference rhf"
		else
			f.puts "set reference uhf"
		end

		# Density fitting - SCF
		if @settings[:density_fitting] == :scf or @settings[:density_fitting] == :both
			# DF setup
			f.puts "set scf_type df"
		else
			f.puts "set df_scf_guess false" unless @settings[:psi4_df_guess]
			f.puts "set scf_type pk"
		end

		# Basis set
		f.puts "basis {"
		# AO basis
		f.puts "   assign #{@settings[:basisset]}"
		@settings.elements_hash(:basisset_elements).each_pair{|element, basis|
			f.puts "   assign #{element.to_s} #{basis}"
		}
		# Auxiliary basis SCF
		if @settings.set?(:auxiliary_basis_scf)
			f.puts "   assign #{@settings[:auxiliary_basis_scf]}-JKFIT df_basis_scf"
		end
		if @settings.set?(:auxiliary_basis_scf_elements)
			@settings.elements_hash(:auxiliary_basis_scf_elements).each_pair{|element, basis|
				f.puts "   assign #{element.to_s} #{basis}-JKFIT df_basis_scf"
			}
		end
		# Auxiliary basis MP2
		if @settings.set?(:auxiliary_basis_mp2)
			f.puts "   assign #{@settings[:auxiliary_basis_mp2]}-RI df_basis_mp2"
			f.puts "   assign #{@settings[:auxiliary_basis_mp2]}-RI df_basis_sapt" if @settings[:method] == :sapt
			f.puts "   assign #{@settings[:auxiliary_basis_mp2]}-RI df_basis_cc" if [:ccsd, :"ccsd(t)"].include?(@settings[:method])
		end
		if @settings.set?(:auxiliary_basis_mp2_elements)
			@settings.elements_hash(:auxiliary_basis_mp2_elements).each_pair{|element, basis|
				f.puts "   assign #{element.to_s} #{basis}-RI df_basis_mp2"
				f.puts "   assign #{element.to_s} #{basis}-RI df_basis_sapt" if @settings[:method] == :sapt
				f.puts "   assign #{element.to_s} #{basis}-RI df_basis_cc" if [:ccsd, :"ccsd(t)"].include?(@settings[:method])
			}
		end
		f.puts "}"

		# Convergence criteria (use default if not specified in cuby)
		f.puts "set e_convergence #{@settings[:scf_convergence]}" if @settings.set?(:scf_convergence)
		f.puts "set d_convergence #{@settings[:density_convergence]}" if @settings.set?(:density_convergence)

		# Threshold for FNO
		f.puts "set occ_tolerance #{@settings[:fno_occ_tolerance]}" if @settings[:psi4_fno]

		# Calculation setup
		f.puts "set guess sad"
		if @settings[:correlation_frozen_core]
			f.puts "set freeze_core true"
		else
			f.puts "set freeze_core false"
		end

		what = "energy"
		what = "gradient" if @what.include?(:gradient)

		# Frozen Natural Orbitals, density fitting in correlation for MP3 and CC methods
		if [:mp3, :ccsd, :"ccsd(t)"].include?(@settings[:method])
			prefix_corr = ""
			prefix_corr = "fno-" if @settings[:psi4_fno]
			if @settings[:density_fitting] == :correlation
				Cuby::error "Psi4 does not allow density-fitted correlation calculation after non-DF SCF."
			elsif @settings[:density_fitting] == :both
				prefix_corr += "df-"
			end
		end

		# Calculation method
		case @settings[:method]
		when :hf
			f.puts "#{what}('scf')"
		when :dft
			require "classes/calculation/dft_options.rb"
			dft_file = File.open(interface_dir + "/dft_options.yaml", "r")
			dft_options = YAML.load(dft_file)
			dft_file.close
			f.puts "#{what}('#{dft_options.functional(@settings)}')"
		when :mp2
			f.puts "set mp2_type conv" unless  @settings[:density_fitting] == :correlation or @settings[:density_fitting] == :both
			f.puts "#{what}('mp2')"
		when :mp3
			f.puts "#{what}('#{prefix_corr}mp3')"
		when :ccsd
			f.puts "#{what}('#{prefix_corr}ccsd')"
		when :"ccsd(t)"
			f.puts "#{what}('#{prefix_corr}ccsd(t)')"
		when :sapt
			# Charge transfer SAPT only
			f.puts "energy('#{@settings[:psi4_sapt_level]}#{@settings[:psi4_sapt_ct] ? '-ct' : ''}')"
		end

		# Print energies as YAML
		# SCF energy
		f.puts "print \":scf_energy:\", get_variable('SCF TOTAL ENERGY')"
		# MP2
		if [:mp2, :mp3, :ccsd, :"ccsd(t)"].include?(@settings[:method])
			f.puts "print \":correlation_mp2:\", get_variable('MP2 CORRELATION ENERGY')"
			f.puts "print \":correlation_mp2_os:\", get_variable('MP2 OPPOSITE-SPIN CORRELATION ENERGY')"
			f.puts "print \":correlation_mp2_ss:\", get_variable('MP2 SAME-SPIN CORRELATION ENERGY')"
		end
		# MP3
		if @settings[:method] == :mp3
			f.puts "print \":correlation_mp3:\", get_variable('MP3 CORRELATION ENERGY')"
		end
		# CCSD
		if [:ccsd, :"ccsd(t)"].include?(@settings[:method])
			f.puts "print \":correlation_ccsd:\", get_variable('CCSD CORRELATION ENERGY')"
			f.puts "print \":correlation_ccsd_os:\", get_variable('CCSD OPPOSITE-SPIN CORRELATION ENERGY')"
			f.puts "print \":correlation_ccsd_ss:\", get_variable('CCSD SAME-SPIN CORRELATION ENERGY')"
			#f.puts "print \"::\", get_variable('')"
		end
		# CCSD(T)
		if @settings[:method] == :"ccsd(t)"
			f.puts "print \":correlation_ccsd(t):\", get_variable('CCSD(T) CORRELATION ENERGY')"
		end

		if @settings[:method] ==:sapt
			f.puts "print \":sapt_ct:\", get_variable('SAPT CT ENERGY')"
			f.puts "print \":sapt_disp:\", get_variable('SAPT DISP ENERGY')"
			f.puts "print \":sapt_elst:\", get_variable('SAPT ELST ENERGY')"
			f.puts "print \":sapt_energy:\", get_variable('SAPT ENERGY')"
			f.puts "print \":sapt_exch:\", get_variable('SAPT EXCH ENERGY')"
			f.puts "print \":sapt_ind:\", get_variable('SAPT IND ENERGY')"
			#--------------------------------------------------
			# f.puts "print \":sapt_:\", get_variable('SAPT SAPT0 ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SAPT2 ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SAPT2+ ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SAPT2+(3) ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SAPT2+(3)(CCD) ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SAPT2+(CCD) ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SAPT2+3 ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SAPT2+3(CCD) ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SCS-DISP ENERGY')"
			# f.puts "print \":sapt_:\", get_variable('SAPT SCS-SAPT0 ENERGY')"
			#-------------------------------------------------- 
		end

		f.close
	end

	def psi4_write_input_geo
		# Write input file using the existing header, only geometry is added
		filename = in_calc_dir("psi4.in")
		f = File.open(filename,"w+")

		# Template from file, substitute geometry
		IO.readlines(in_calc_dir("psi4.header")).each{|line|
			if line =~ /%geometry%/
				# Write geometry')"

				if @settings[:psi4_geometry_fragments]
					# Fragments written separately
					psi4_write_geo_fragments(f)
				else
					psi4_write_geo_block(f, @geometry, @settings[:charge], @settings[:multiplicity])
				end
			else
				# Write the original line
				f.puts line
			end
		}

		f.close
	end

	def psi4_write_geo_fragments(file)
		# Assumes there are two blocks with molecule definition
		[:molecule_a, :molecule_b].each{|molname|
			file.puts "   --" unless molname == :molecule_a
			geo = @geometry.geometry_from_selection(@settings[molname, :selection])
			charge = @settings[molname, :charge]
			multiplicity = @settings[molname, :multiplicity]
			psi4_write_geo_block(file, geo, charge, multiplicity)
		}
	end

	def psi4_write_geo_block(file, geo, charge, multiplicity)
		file.puts "   #{charge} #{multiplicity}"
		geo.each_index{|i|
			atom = geo[i]
			element = atom.element.to_s
			element = "@" + element if atom.properties[:ghost]
			file.printf("%5s%25.15f%25.15f%25.15f\n", element, atom.x, atom.y, atom.z)
		}
	end

	def psi4_run
		# Build input with the current geometry
		psi4_write_input_geo

		# Run psi4
		command = "cd #{calc_dir};"
		command << "export PATH=#{@settings[:psi4_bin_dir]}:$PATH;"
		command << "export PSIDATADIR=#{@settings[:psi4_data_dir]};" unless @settings[:psi4_data_dir] == ""
		command << "export LC_ALL=C;"
		command << "export PSI_SCRATCH=`pwd`;"
		command << "psi4 -n #{@settings[:parallel]} -i psi4.in -o psi4.out > output.yaml 2> psi4.err;"
		unless system(command)
			# Psi4 failed - read common errors from the standard output
			parser = OutputParser.new(in_calc_dir("psi4.out"), @name)
			parser.add_error(/Unable to load (.*-jkfit)\.gbs from the default Psi4 basis set library/, "Auxiliary basis '%1' not found in Psi4 library")
			parser.add_error(/Unable to load (.*-mp2fit)\.gbs from the default Psi4 basis set library/, "Auxiliary basis '%1' not found in Psi4 library")
			parser.add_error(/Unable to load (.*)\.gbs from the default Psi4 basis set library/, "Basis set '%1' not found in Psi4 library")
			parser.add_error(/PsiException: Keyword DF_BASIS_MP2 is required/, "Auxiliary basis for DF-MP2 can not be determined automatically")
			parser.add_error(/error: not enough memory for ccsd/, "Not enough memory for CCSD")
			#parser.add_error(//, "")
			parser.execute

			# Psi4 failed - read common errors from the stderr output
			parser = OutputParser.new(in_calc_dir("psi4.err"), @name)
			parser.add_error(/RuntimeError: Not enough memory/, "Not enough memory")
			#parser.add_error(//, "")
			parser.execute

			#Generic error
			Cuby::error "PSI4 returned nonzero exit code (calculation #{@name})"
		end
	end

	def psi4_read # => Results
		results = Results.new
		
		# Energies are read from a separate YAML file
		f = File.open(in_calc_dir("output.yaml"))
		YAML.load(f).each_pair{|key, value|
			results.energy_components[key] = value * HARTREE2KCAL
		}

		results.energy = results.energy_components[:scf_energy]

		case @settings[:method]
		when :mp2
			results.energy += results.energy_components[:correlation_mp2]
		when :mp3
			results.energy += results.energy_components[:correlation_mp3]
		when :ccsd
			results.energy += results.energy_components[:correlation_ccsd]
		when :"ccsd(t)"
			results.energy += results.energy_components[:"correlation_ccsd(t)"]
		when :sapt
			results.energy = results.energy_components[:sapt_energy]
		end

		# Gradient is read from the output file
		if @what.include?(:gradient)
			parser = OutputParser.new(in_calc_dir("psi4.out"), @name)
			results.gradient = Gradient.new

			parser.add_block(
				:geometry,
				/^ +Geometry \(in Angstrom\)/, 
				/^ +\S+ +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+)/, 
				/Nuclear repulsion/,
				{:get => :match_array, :type => :float}
			)
			parser.add_block(
				:gradient,
				/^ +-Total [G,g]radient/, 
				/^ +[0-9]+ +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+)/, 
				/^$/,
				{:get => :match_array, :type => :float}
			)

			geo_new = Geometry.new
			gi = 0
			parser.execute{|name, result, count|
				case name
				when :geometry
					x,y,z = result
					geo_new << Atom.new(@geometry[gi].element,x,y,z)
					gi += 1
				when :gradient
					x,y,z = result
					results.gradient << Coordinate[x,y,z] * HARTREE2KCAL / BOHR2ANGSTROM
				end
			}
		end

		return results
	end

end
