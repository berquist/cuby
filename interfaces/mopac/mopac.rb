################################################################################
#
# MOPAC interface
#
# Author: Jan Rezac
# Date created: 2013-04-09
# License: Cuby4 license
# Description: Interface for external calculations in MOPAC
# Status: Works
#
################################################################################

#===============================================================================
# Interface to the MOPAC 2012 (works also with version 2009)
# Free for academic use, available after registration
# http://www.openmopac.net
#===============================================================================

require "classes/tools/output_parser.rb"

module InterfaceMopac
	#=======================================================================
	# Interface header
	#=======================================================================
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "All cuby3 functionality implemented, examples for all features"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :atomic_charges]
	MODIFIER = false
	DIRECTORY = "MOPAC"
	# Methods provided by the interface:
	METHODS = {
		:"mndo"		=> [:energy, :gradient, :atomic_charges],
		:"mndo-d"	=> [:energy, :gradient, :atomic_charges],
		:"am1" 		=> [:energy, :gradient, :atomic_charges],
		:"rm1" 		=> [:energy, :gradient, :atomic_charges],
		:"pm3" 		=> [:energy, :gradient, :atomic_charges],
		:"pm6" 		=> [:energy, :gradient, :atomic_charges],
		:"pm7" 		=> [:energy, :gradient, :atomic_charges]
	}
	SOLVENTS = [:cosmo]
	ATOMIC_CHARGES = [:mulliken]
	#=======================================================================


	# Interface constants (lookup tables for construction of the input)
	MOPAC_MULTIPLICITY = {
		1 => "SINGLET",
		2 => "DOUBLET",
		3 => "TRIPLET",
		4 => "QUARTET",
		5 => "QUINTET"
	}
	MOPAC_METHODS = {
		:"mndo" 	=> "MNDO",
		:"mndo-d" 	=> "MNDOD",
		:"am1" 		=> "AM1",
		:"rm1" 		=> "RM1",
		:"pm3" 		=> "PM3",
		:"pm6" 		=> "PM6",
		:"pm7" 		=> "PM7"
	}
	MOPAC_CORRECTIONS = {
		:none =>	"",
		:"d2" =>	"-D2",
		:"h2" =>	"-H2",
		:"dh2" =>	"-DH2",
		:"dh+" =>	"-DH+",
		:"d3" =>	"-D3",
		:"h4" =>	"-H4",
		:"d3h4" =>	"-D3H4",
		:"d3h4x" =>	"-D3H4X"
	}


	def prepare_interface
		# Compatibility checks
		if @settings[:mopac_version] == :mopac7
			if @geometry.size <= 3
				# For small molecules, mopac7 can use only internal coordinates
				Cuby::error("Mopac7 interface supports only molecules with 4 or more atoms")
			end
		end

		# Things that have to be done even when old input is reused:

		# Setpi - manual definition of double bonds for mozyme
		if @settings[:mopac_setpi].size > 0
			@setpi_list = [] # Array of the lines added after geometry
			@settings[:mopac_setpi].each{|rec|
				rec.gsub!(";",",") # For backward compatibility when selections were not used
				list = @geometry.atomlist_from_selection(rec)
				Cuby::error "SETPI: two atoms should be selected (now #{list.size})" unless list.size == 2
				@setpi_list << list.map{|x| (x+1).to_s}.join(" ")
			}
		end

		# Manual assignment of atomic charges
		@charge_selections = {}
		@settings[:mopac_setcharge].each_pair{|sel, value|
			list = @geometry.atomlist_from_selection(sel.to_s)
			list.each{|i| @charge_selections[i] = value}
		}

		# Prepare calculation directory
		if calc_dir_mkdir("mopac.header", "mopac.in.out") == :old
			calc_using_input # Old input should be reused
			return
		end

		# Write the keyword part of the input
		mopac_write_keywords # This writes file mopac.header

		# Create complete input by combining mopac.header and the current geometry
		mopac_write_input_geo

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		mopac_run unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/mopac.in.out")
		return mopac_read
	end

	def cleanup_interface
		calc_dir_delete if @settings[:job_cleanup]
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def mopac_write_keywords
		#: Write the header of the input where the calculation is set up
		#: geometry is added later, so that the header can be modified
		#: by the user after the job is prepared

		filename = in_calc_dir("mopac.header")
		f = File.open(filename,"w+")

		# Method
		f.print MOPAC_METHODS[@settings[:method]]

		# Corrections
		f.print MOPAC_CORRECTIONS[@settings[:mopac_corrections]]

		# SCF
		f.print " 1SCF DENOUT"
		f.print " CHARGE=#{@settings[:charge]}"

		# Spin
		f.printf " #{MOPAC_MULTIPLICITY[@settings[:multiplicity]]}"
		if @settings[:spin_restricted]
			f.printf " RHF" if @settings[:mopac_version] == :mopac2012
		else
			f.printf " UHF"
		end

		# Mozyme
		f.print " MOZYME" if @settings[:mopac_mozyme]

		# Setpi - manual definition of double bonds for mozyme
		if @settings[:mopac_setpi].size > 0
			# Write to input
			f.print " SETPI"
		end

		# Precise - accurate gradient
		f.print " PRECISE" if @settings[:mopac_precise]

		# COSMO
		if @settings[:solvent_model] == :cosmo
			f.print " EPS=#{@settings[:solvent_epsilon]}"
		end

		# Calculate gradient or charges if needed
		#f.print " GRADIENTS DCART" if @what.include?(:gradient)
		if @what.include?(:gradient)
			f.print " GRADIENTS AUX(MOS=-99999,XW,XS,XP,PRECISION=6)"
		end
		f.print " MULLIK" if @what.include?(:atomic_charges)

		# Parallelization
		if @settings[:mopac_parallel] && @settings[:mopac_version] == :mopac2012
			f.print " THREADS=#{@settings[:parallel]}"
		end

		# Correction for peptide bonds (addition for RM1, PM6, and PM7)
		if @settings[:mopac_peptide_bond_fix] && [:pm6, :pm7, :rm1].include?(@settings[:method])
			f.print " MMOK"
		else
			f.print " NOMM"
		end

		# Print the corrections included in PM7
		if @settings[:method] == :pm7
			f.print " DISP"
		end


		# Extra keywords
		f.print " " + @settings[:mopac_keywords] if @settings[:mopac_keywords] != ""

		f.close
	end

	def mopac_write_input_geo
		#: Write final input for mopac, using the previously prepared header
		#: and a current geometry

		# Read mopac.header
		header = IO.readlines(in_calc_dir("mopac.header"))

		# Complete input:
		filename = in_calc_dir("mopac.in")
		f = File.open(filename,"w+")

		# Print header from the header file, adding option to restart from previous density
		restart = ""
		if FileTest.exist?(in_calc_dir("mopac.in.den")) && @settings[:start_from_previous]
			restart =  " OLDENS" 
		end
		f.puts header[0].strip + restart

		# Separator
		f.print "\n\n"
		case @settings[:mopac_version]
		when :mopac2012
			@geometry.each_index{|i|
				atom = @geometry[i]
				ch = ""
				if @charge_selections[i]
					ch = "(#{@charge_selections[i]})"
				end
				f.printf("  %2s%-4s    %15.10f %15.10f %15.10f\n",atom.element,ch,atom.x,atom.y,atom.z)
			}
		when :mopac7
			@geometry.each_index{|i|
				atom = @geometry[i]
				f.printf("  %2s   %15.10f 0 %15.10f 0 %15.10f 0\n",atom.element,atom.x,atom.y,atom.z)
			}
		end
		f.puts

		# Setpi - manual definition of double bonds for mozyme
		# Write the list to the input
		if @setpi_list
			@setpi_list.each{|line| f.puts line}
			f.puts
		end

		f.close
	end

	def mopac_run
		#: Set up environment for MOPAC and run it

		# Build input with the current geometry
		mopac_write_input_geo

		# Prepare the command for running mopac
		command = "cd #{calc_dir}; "
		case @settings[:mopac_version]
		when :mopac2012
			command << "export MOPAC_LICENSE=#{File.dirname(@settings[:mopac_exe])};"
			command << "export LD_LIBRARY_PATH=#{File.dirname(@settings[:mopac_exe])}:$LD_LIBRARY_PATH;"
			command << "#{@settings[:mopac_exe]} mopac.in 2> mopac.err;"
		when :mopac7
			command << "export FOR005=mopac.in; "
			command << "export FOR006=mopac.out2; "
			command << "export FOR009=mopac.res; "
			command << "export FOR010=mopac.den; "
			command << "export FOR011=mopac.log; "
			command << "export FOR012=mopac.arc; "
			command << "#{@settings[:mopac_exe]} 2> mopac.err > mopac.in.out"
		end

		# Run mopac
		unless system(command)
			Cuby::error "Mopac returned nonzero exit code (calculation #{@name})"
		end
	end

	def mopac_read #u=> Results
		#: Reads energy, gradient and dipole from MOPAC output
	
		results = Results.new
		parser = OutputParser.new(in_calc_dir("mopac.in.out"), @name)

		# Energy
		parser.add_pattern(:energy, /^ +FINAL HEAT OF FORMATION = +([-]*[0-9]*\.[0-9]+)/, :get => :submatch)

		if @settings[:method] == :pm7
			parser.add_pattern(:energy_d,  /^ +DISPERSION ENERGY += +([-]*[0-9]*\.[0-9]+)/, :get => :submatch, :required => false)
			parser.add_pattern(:energy_hb, /^ +HYDROGEN-BOND ENERGY += +([-]*[0-9]*\.[0-9]+)/, :get => :submatch, :required => false)
		end

		# Dipole
		parser.add_pattern(:dipole,
			/^ SUM\s+([-]*[0-9]*\.[0-9]+)\s+([-]*[0-9]*\.[0-9]+)\s+([-]*[0-9]*\.[0-9]+)\s+([-]*[0-9]*\.[0-9]+)/,
			:get => :match_array, :type => :float, :required => false
		)

		# Errors
		parser.add_error(/CHARGE SPECIFIED FOR SYSTEM .* IS INCORRECT/, "Mopac thinks the charge is wrong")

		# Add charges reader if needed
		if what.include?(:atomic_charges)
			results.atomic_charges = AtomicCharges.new(:mulliken)
			# The default charges printed by MOPAC
			#--------------------------------------------------
			# parser.add_block(
			# 	:charges, 
			# 	/NET ATOMIC CHARGES/, 
			# 	/^\s+[0-9]+\s+[A-Z,a-z]+\s+([-]*[0-9]*\.[0-9]+)\s+/,
			# 	/^ DIPOLE/,
			# 	{:get => :submatch, :type => :float}
			# )
			#-------------------------------------------------- 

			# Mulliken charges
			parser.add_block(
				:charges, 
				/MULLIKEN POPULATIONS AND CHARGES/, 
				/^\s+[0-9]+\s+[A-Z,a-z]+\s+[^ ]+\s+([-]*[0-9]*\.[0-9]+)\s+/,
				/^ TOTAL CPU TIME/,
				{:get => :submatch, :type => :float}
			)
		end

		# Run parser
		parser.execute{|name, result, count|
			case name
			when :charges
				results.atomic_charges << result
			end
		}

		# Energy
		results.energy = parser[:energy]

		# PM7 corrections
		if @settings[:method] == :pm7
			results.energy_components[:dispersion] = parser[:energy_d] if parser[:energy_d]
			if parser[:energy_hb]
				results.energy_components[:h_bonds] = parser[:energy_hb] 
			else
				results.energy_components[:h_bonds] = 0.0
			end

		end

		# Dipole
		if parser[:dipole]
			results.multipoles[:dipole] = Multipoles::Dipole.new(
				parser[:dipole][0] * DEBYE2CUBY,
				parser[:dipole][1] * DEBYE2CUBY,
				parser[:dipole][2] * DEBYE2CUBY
			)
		end

		# Gredient is read from the .aux file
		if what.include?(:gradient)
			mopac_aux = MopacAux.new(in_calc_dir("mopac.in.aux"))
			results.gradient = mopac_aux.gradient
		end

		return results
	end

	class MopacAux
		attr_reader :gradient

		def initialize(filename)
			# Read and parse the aux file

			f = File.open(filename, "r")

			# Block readers
			read_gradient = false
			gradient = []

			# Iterate through file
			while line = f.gets
				# Remove spaces around string
				line.strip!
				# Find blocks
				if line =~ /^GRADIENTS:KCAL\/MOL\/ANGSTROM/
					read_gradient = true
					next
				end
				# Terminate blocks
				if line[0..0] =~ /[A-Z,a-z]/
					read_gradient = false
					next
				end
				# Gradient block
				if read_gradient
					line.gsub("-"," -").split.each{|x| gradient << x.to_f}
				end
			end
			f.close

			# Build the results
			@gradient = Gradient.from_vector(gradient)
		end
	end
end
