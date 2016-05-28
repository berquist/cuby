################################################################################
#
# ORCA interface
#
# Author: Jan Rezac
# Date created: 2013-04-29
# License: Cuby4 license
# Description: Interface for external calculations using ORCA
# Status: Works
#
################################################################################

#===============================================================================
# Interface to ORCA
# http://www.cec.mpg.de/forum/portal.php
#===============================================================================

# Important notes:
#   * Geometry can not be readed from an external file, ghost atoms can not be
#     defined in such case.

require "classes/tools/output_parser.rb"

module InterfaceOrca
	#=======================================================================
	# Interface header
	#=======================================================================
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "Works"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :ghost_atoms, :point_charges, :point_charges_gradient]
	MODIFIER = false
	DIRECTORY = "ORCA"
	# Methods provided by the interface:
	METHODS = {
		:"hf"		=> [:energy, :gradient, :point_charges, :point_charges_gradient],
		:"dft"		=> [:energy, :gradient, :point_charges, :point_charges_gradient],
		:"mp2"		=> [:energy, :gradient],
		:"mp3"		=> [:energy],
		:"ccsd"		=> [:energy],
		:"ccsd(t)"	=> [:energy],
		:"qcisd"	=> [:energy],
		:"qcisd(t)"	=> [:energy]
	}
	SOLVENTS = [:cosmo]
	# Unsupported keywords: options common in similar interfaces that can't
	# be used
	UNSUPPORTED_KEYWORDS = [
		:use_symmetry # Symmetry is not supported meaningfully in orca
	]
	#=======================================================================

	@@orca_basis_sets = nil

	def prepare_interface
		# Load basis set list on the first use
		unless @@orca_basis_sets
			File.open(interface_dir + "/basissets.yaml"){|f|
				@@orca_basis_sets = YAML.load(f)
			}
		end

		# Prepare calculation directory
		if calc_dir_mkdir("orca.header", "orca.out") == :old
			calc_using_input # Old input should be reused
			return
		end

		# Write the keyword part of the input
		orca_write_keywords # This writes file orca.header

		# Create complete input by combining orca.header and the current geometry
		orca_write_input_geo

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		orca_run unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/orca.out")
		return orca_read
	end

	def cleanup_interface
		if @settings[:job_cleanup]
			# Delete everything
			calc_dir_delete 
		elsif @settings[:delete_large_files]
			# Delete only large files (and some unnecessary files as well)
			system("cd #{calc_dir}; rm -f orca.gbw")
		end
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def orca_write_keywords
		# Write header of the input file
		filename = in_calc_dir("orca.header")
		f = File.open(filename,"w+")

		# Density fitting
		dft_ri = @settings[:density_fitting] == :scf || @settings[:density_fitting] == :both
		correlation_ri = @settings[:density_fitting] == :correlation || @settings[:density_fitting] == :both

		# Memory
		mem = @settings[:mem] * @settings[:mem_core] / 100
		f.puts "%MaxCore #{mem}"

		# Paralellization
		if @settings[:parallel] > 1
			f.puts "%pal"
		        f.puts " nprocs #{@settings[:parallel]}"
			f.puts "end"
		end

		# Method
		method = ""

		# Method prefix
		if @settings[:orca_dlpno] # DLPNO
			# Checks
			Cuby::error "DLPNO available only with CCSD and CCSD(T)" unless [:ccsd, :"ccsd(t)"].include?(@settings[:method])
			Cuby::error "DLPNO: keyword auxiliary_basis_mp2 must be set" unless @settings.set?(:auxiliary_basis_mp2)


			method += "DLPNO-"
			rimp2_basis = "/C"
		end

		# Method selection
		case @settings[:method]
		when :hf
			method += "HF"
		when :dft
			method += "DFT"
			# Load DFT setup
			require "classes/calculation/dft_setup.rb"
			dft_setup = DftSetup.new(interface_dir + "/dft_functionals.yaml", interface_dir + "/dft_grids.yaml")

			# Determin RI mode
			riscf_basis = nil
			if dft_ri
				case dft_setup.functional_type(@settings)
				when :gga
					method += " RI"
					riscf_basis = "/J"
				when :meta_gga
					method += " RI"
					riscf_basis = "/J"
				when :hybrid
					method += " RIJK"
					riscf_basis = "/JK"
				when :range_separated
					Cuby::error "Orca interface: RI not supported for range-separated functionals"
				else
					Cuby::error "Cuby does not implement RI settings for #{dft_setup.functional_type(@settings)} functionals"
				end
			end
			
			# Functional
			method += " #{dft_setup.functional(@settings)}"
		when :mp2
			method += "MP2"
			rimp2_basis = "/C" if correlation_ri
		when :mp3
			method += "MP3"
		when :ccsd
			method += "CCSD"
		when :"ccsd(t)"
			method += "CCSD(T)"
		when :qcisd
			method += "QCISD"
		when :"qcisd(t)"
			method += "QCISD(T)"
		end



		# RHF/UHF
		if @settings[:spin_restricted]
			method += " RHF"
		else
			method += " UHF"
		end

		# SCF convergence
		method += " SCFCONV#{@settings[:scf_convergence]}"

		# Restart
		method += " NoAutoStart" if !@settings[:start_from_previous]

		f.puts "# Orca calculation for cuby"
		f.puts "! #{method}"

		# Extra keywords, "!" added if needed
		if @settings.set?(:orca_keywords)
			line = @settings[:orca_keywords]
			if line =~ /^!/
				f.puts line
			else
				f.puts "! " + line
			end
		end

		# Point charges
		if @point_charges
			f.puts '% pointcharges "pointcharges.pc"'
		end

		# Basis set
		check_basis(@settings[:basisset])
		basis = @settings[:basisset]

		# Auxiliary basis
		if @settings.set?(:auxiliary_basis_scf) && riscf_basis
			check_auxiliary_basis(@settings[:auxiliary_basis_scf] + riscf_basis)
			basis += " #{@settings[:auxiliary_basis_scf]}#{riscf_basis}"
		end
		if @settings.set?(:auxiliary_basis_mp2) && rimp2_basis
			check_auxiliary_basis(@settings[:auxiliary_basis_mp2] + rimp2_basis)
			basis += " #{@settings[:auxiliary_basis_mp2]}#{rimp2_basis}"
		end

		case @settings[:orca_ecp]
		when :auto
			f.puts "! ECP{#{basis.gsub(' ',',')}} #{basis}"
		when :none
			f.puts "! #{basis}"
		end

		# Basis block
		f.puts "%basis"
		@settings.elements_hash(:basisset_elements).each_pair{|element, basis|
			check_basis(basis)
			f.puts "   NewGTO #{PeriodicTable.proton_number(element)}"
			f.puts "     \"#{basis}\""
			f.puts "   end"
			if riscf_basis
				unless @settings.elements_hash(:auxiliary_basis_scf_elements)[element]
					Cuby::error "Non-default basis is set for element #{element}, but SCF auxiliary basis is not defined for it"
				end
			end
			if rimp2_basis
				unless @settings.elements_hash(:auxiliary_basis_mp2_elements)[element]
					Cuby::error "Non-default basis is set for element #{element}, but MP2 auxiliary basis is not defined for it"
				end
			end

		}
		# Auxiliary basis per element
		if dft_ri
			@settings.elements_hash(:auxiliary_basis_scf_elements).each_pair{|element, basis|
				check_auxiliary_basis("#{basis}#{riscf_basis}")
				f.puts "   NewAuxGTO #{PeriodicTable.proton_number(element)}"
				f.puts "     \"#{basis}#{riscf_basis}\""
				f.puts "   end"
			}
		end
		if rimp2_basis
			@settings.elements_hash(:auxiliary_basis_mp2_elements).each_pair{|element, basis|
				check_auxiliary_basis("#{basis}#{rimp2_basis}")
				f.puts "   NewAuxGTO #{PeriodicTable.proton_number(element)}"
				f.puts "     \"#{basis}#{rimp2_basis}\""
				f.puts "   end"
			}
		end
		f.puts "end"


		# Method setup % method
		f.puts "%method"
		f.puts "   RunTyp Gradient" if @what.include?(:gradient) #!# not tested

		# DFT
		if @settings[:method] == :dft
			# Functional setup (the following line) was moved to the !-line at the beginning of the input
			# f.puts "   Functional #{dft_setup.functional(@settings)}"

			if @settings[:dft_grid] == :custom
				f.puts "   Grid #{@settings[:dft_grid_custom]}"
			else
				f.puts "   Grid #{dft_setup.grid(@settings)[0]}"
				f.puts "   FinalGrid #{dft_setup.grid(@settings)[1]}"
			end
		end

		# MP2
		if @settings[:method] == :mp2
			f.puts @settings[:correlation_frozen_core] ? "   FrozenCore FC_ELECTRONS" : "   FrozenCore FC_NONE"
		end

		# End of %method
		f.puts "end"

		# MP2 setup
		if @settings[:method] == :mp2
			# MP2 block
			f.puts "%mp2"
			# RI on/off
			f.puts correlation_ri ? "   RI on" : "   RI off"
			# Print more results
			f.puts "   PrintLevel 3"
			# SCS-MP2
			if @settings[:spin_component_scaling]
				f.puts "   DoSCS true"
				f.puts "   Ps #{@settings[:scs_os]}"
				f.puts "   Pt #{@settings[:scs_ss]}"
			end
			# End of MP2 block
			f.puts "end"
		end

		# COSMO
		if @settings[:solvent_model] == :cosmo
			f.puts "%cosmo"
			f.puts "   epsilon #{@settings[:solvent_epsilon]}"
			f.puts "end"
		end

		# Extra input - without any processing
		if @settings.set?(:orca_extra_input)
			f.puts @settings[:orca_extra_input]
		end

		# Coordinates start
		f.puts "* xyz #{@settings[:charge]} #{@settings[:multiplicity]}"
		f.close
	end

	def orca_write_input_geo
		# Write input file using the existing header, only geometry is added
		filename = in_calc_dir("orca.in")
		f = File.open(filename,"w+")

		# Header from file
		header = IO.readlines(in_calc_dir("orca.header"))
		f.puts header

		# Geometry
		@geometry.each_index{|i|
			atom = @geometry[i]
		 	ghost = atom.properties[:ghost] ? ":" : " "
		 	f.printf(" %3s %s %15.10f %15.10f %15.10f\n", atom.element.to_s, ghost, atom.x, atom.y, atom.z)
		}
		f.puts "*"

		f.close
	end

	def orca_run
		# Build input with the current geometry
		orca_write_input_geo
		# Point charges
		orca_write_point_charges if @point_charges

		# Run orca
		command = "export PATH=#{@settings[:orca_bin_dir]}:$PATH;"
		command << @settings[:orca_mpi_setup]
		command << "cd #{calc_dir};"
		command << "#{@settings[:orca_bin_dir]}/orca orca.in > orca.out 2> orca.err;"
		unless system(command)
			Cuby::error "ORCA returned nonzero exit code (calculation #{@name})"
		end
	end

	def orca_read # => Results
		results = Results.new
		parser = OutputParser.new(in_calc_dir("orca.out"), @name)

		# SCF energy
		parser.add_pattern(:e_scf, /^Total Energy       : *([^ ]*)/, :get => :submatch)

		# Method-specific
		case @settings[:method]
		when :mp2
			if @settings[:density_fitting] == :correlation || @settings[:density_fitting] == :both
				parser.add_pattern(:corr_mp2, /^ RI-MP2 CORRELATION ENERGY/)
			else
				parser.add_pattern(:corr_mp2, /^ MP2 CORRELATION ENERGY/)
			end
		when :mp3
			parser.add_pattern(:corr_mp2, /EC\(MP2\)= *([^ ]*)/, :get => :submatch)
			parser.add_pattern(:corr_mp3, /EC\(MP3\)= *([^ ]*)/, :get => :submatch)
		when :"ccsd"
			parser.add_pattern(:"e_ccsd", /^FINAL SINGLE POINT ENERGY/)
		when :"ccsd(t)"
			parser.add_pattern(:"e_ccsd", /^E\(CCSD\)/)
			parser.add_pattern(:"e_ccsd(t)", /^E\(CCSD\(T\)\)/)
		when :qcisd
			parser.add_pattern(:"e_qcisd", /^FINAL SINGLE POINT ENERGY/)
		when :"qcisd(t)"
			parser.add_pattern(:"e_qcisd", /^E\(QCISD\)/)
			parser.add_pattern(:"e_qcisd(t)", /^E\(QCISD\(T\)\)/)
		end

		# Gradient
		if @what.include?(:gradient)
			results.gradient = Gradient.new
			if @settings[:method] == :mp2
				parser.add_block(
					:gradient,
					/The final MP2 gradient/, 
					/^[^:]+: +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+)/, 
					/NORM/,
					{:get => :match_array, :type => :float}
				)
			else
				parser.add_block(
					:gradient,
					/CARTESIAN GRADIENT/, 
					/^[^:]+: +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+) +([-]*[0-9]*\.[0-9]+)/, 
					/Norm/,
					{:get => :match_array, :type => :float}
				)
			end
		end

		# DLPNO
		if @settings[:orca_dlpno]
			parser.add_pattern(:corr_sl_pno_mp2, /SL-PNO-MP2 correlation energy \(all kept pairs\)[. ]*([-]*[0-9]*\.[0-9]+)/, :get => :submatch)
			parser.add_pattern(:corr_dlpno_init, /Initial total correlation energy[. ]*([-]*[0-9]*\.[0-9]+)/, :get => :submatch)
		end

		# Run parser
		parser.execute{|name, result, count|
			case name
			when :gradient
				x,y,z = result
				results.gradient << Coordinate[x,y,z] * HARTREE2KCAL / BOHR2ANGSTROM
			end
		}

		# SCF energy
		results.energy = parser[:e_scf] * HARTREE2KCAL

		# Correlation
		case @settings[:method]
		when :mp2
			results.energy_components[:scf_energy] = results.energy
			if @settings[:spin_component_scaling]
				results.energy_components[:correlation_scs_mp2] = parser[:corr_mp2] * HARTREE2KCAL
				results.energy += results.energy_components[:correlation_scs_mp2]
			else
				results.energy_components[:correlation_mp2] = parser[:corr_mp2] * HARTREE2KCAL
				results.energy += results.energy_components[:correlation_mp2]
			end
		when :mp3
			results.energy_components[:scf_energy] = results.energy
			results.energy_components[:correlation_mp2] = parser[:corr_mp2] * HARTREE2KCAL
			results.energy_components[:correlation_mp3] = parser[:corr_mp3] * HARTREE2KCAL
			results.energy += results.energy_components[:correlation_mp3]
		when :"ccsd"
			results.energy_components[:scf_energy] = results.energy
			results.energy_components[:"correlation_ccsd"] = parser[:"e_ccsd"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			results.energy = parser[:"e_ccsd"] * HARTREE2KCAL
		when :"ccsd(t)"
			results.energy_components[:scf_energy] = results.energy
			results.energy_components[:"correlation_ccsd"] = parser[:"e_ccsd"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			results.energy_components[:"correlation_ccsd(t)"] = parser[:"e_ccsd(t)"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			results.energy = parser[:"e_ccsd(t)"] * HARTREE2KCAL
		when :"qcisd"
			results.energy_components[:scf_energy] = results.energy
			results.energy_components[:"correlation_qcisd"] = parser[:"e_qcisd"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			results.energy = parser[:"e_qcisd"] * HARTREE2KCAL
		when :"qcisd(t)"
			results.energy_components[:scf_energy] = results.energy
			results.energy_components[:"correlation_qcisd"] = parser[:"e_qcisd"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			results.energy_components[:"correlation_qcisd(t)"] = parser[:"e_qcisd(t)"] * HARTREE2KCAL - results.energy_components[:scf_energy]
			results.energy = parser[:"e_qcisd(t)"] * HARTREE2KCAL
		end

		if @settings[:orca_dlpno]
			results.energy_components[:correlation_pl_pno_mp2] = parser[:corr_sl_pno_mp2] * HARTREE2KCAL
			results.energy_components[:correlation_dlpno_init] = parser[:corr_dlpno_init] * HARTREE2KCAL
		end

		# Gradient on point charges
		if @what.include?(:gradient) && @point_charges
		        @point_charges.gradient = orca_read_point_charges_gradient
			results.point_charges_gradient = @point_charges.gradient
		end

		return results
	end

	def check_auxiliary_basis(basis_name)
		unless @@orca_basis_sets["fitting"].include?(basis_name.downcase)
			Cuby::recommendation "The selected auxiliay basis is not on the list of basis sets induded in ORCA.\nThis list might be incomplete but check the input before continuing."
		end
	end

	def check_basis(basis_name)
		unless @@orca_basis_sets["ao"].include?(basis_name.downcase)
			Cuby::recommendation "The selected basis set is not on the list of basis sets induded in ORCA.\nThis list might be incomplete but check the input before continuing."
		end
	end

	def orca_write_point_charges
		# Write file with point charges in the calculation directory
		filename = in_calc_dir("pointcharges.pc")
		f = File.open(filename,"w+")
		f.puts @point_charges.size
		@point_charges.each { |pcharge|
			f.printf("%15.8f%15.8f%15.8f%15.8f\n", pcharge.charge, pcharge.x, pcharge.y, pcharge.z)
		}
		f.close
	end

	def orca_read_point_charges_gradient
		# Read gradient on point charges
		gradient = Gradient.new
		filename = in_calc_dir("orca.pcgrad")
		f = File.open(filename,"r")
		ncharges = f.gets.to_i
		unless ncharges == @point_charges.size
			Cuby::error "Number of point charge gradients read does not match input"
		end
		ncharges.times{|i|
			g = f.gets.strip.gsub("D","e").split.map{|x| x.to_f}
			gradient << Coordinate[g[0],g[1],g[2]] * HARTREE2KCAL / BOHR2ANGSTROM
		}
		f.close

		return gradient
	end
end
