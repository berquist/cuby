################################################################################
#
# AMBER interface
#
# Author: Jan Rezac
# Date created: 2012-08-23
# License: Cuby4 license
# Description: Interface for external calculations in AMBER
# Status: Works
#
################################################################################

#===============================================================================
# Interface to the AMBER molecular mechanics package
# Commercial software
# http://ambermd.org/
#===============================================================================

# Possible improvements
# * Option to run the calculations using NAB
# * Automated charge generation using Antechamber (am1-bcc charges first)

require "erb" # ERB templating system, used for construction of the input
require "classes/tools/output_parser.rb"

module InterfaceAmber
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "All cuby3 functionality implemented, examples for all features"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :atomic_charges, :point_charges, :point_charges_gradient]
	MODIFIER = false
	DIRECTORY = "AMBER"
	# Methods provided by the interface:
	METHODS = {
		:forcefield		=> [:energy, :gradient, :atomic_charges, :point_charges, :point_charges_gradient]
	}
	SOLVENTS = [:gbm, :igb7]
	ATOMIC_CHARGES = [:forcefield]
	#=======================================================================


	def prepare_interface
		# Check config
		amber_check_config

		# Check geometry - it must contain the PDB data
		unless @geometry[0].properties[:pdb_atom_name]
			Cuby::error("Amber calculations are possible only on geometry from a PDB file")
		end

		# Prepare calculation directory
		if calc_dir_mkdir("sander.in", "sander.out") == :old
			calc_using_input
			return
		end

		# Write geometry
		GeoFile::AmberCrd.write(@geometry, :file => calc_dir+"/crd")

		# Write input
		if @settings.set?(:amber_input_file)
			FileUtils.cp(@settings[:amber_input_file], calc_dir+"/sander.in")
		else
			amber_write_input(calc_dir+"/sander.in")
		end

		# Generate or copy top file
		if @settings.set?(:amber_top_file)
			FileUtils.cp(@settings[:amber_top_file], calc_dir+"/prmtop")
		else
			amber_generate_prmtop
		end

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		amber_run_calculation unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/sander.out")
		results = amber_read_results
		results.atomic_charges = amber_read_charges if @what.include?(:atomic_charges)
		if @what.include?(:point_charges_gradient)
			results.point_charges_gradient = @point_charges.gradient
		end
		return results
	end

	def cleanup_interface
		calc_dir_delete if @settings[:job_cleanup]
	end

	#=======================================================================
	# Extra methods - MM specific
	#=======================================================================
	
	def get_atomic_charges
		if ! @atomic_charges
			@atomic_charges = amber_read_charges
		end
		return @atomic_charges
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def amber_write_input(filename)
		# Load the template and set required variables
		case @settings[:solvent_model]
		when :gbm
			gbmtype = 1
			template = ERB.new(IO.read(interface_dir + "/templates/sander.in_GBM.erb"))
		when :igb7
			gbmtype = 7
			template = ERB.new(IO.read(interface_dir + "/templates/sander.in_GBM.erb"))
		when :none
			template = ERB.new(IO.read(interface_dir + "/templates/sander.in_vacuum.erb"))
		else
			Cuby::error "Solvent model not supported by AMBER"
		end

		file = File.open(filename,"w+")
		file.print(template.result(binding))
		file.close
	end

	def amber_run_calculation
		# If point charges are used, and atomic charges are not available, get atomic charges
		if @point_charges && ! @atomic_charges
			@atomic_charges = amber_read_charges
		end

		# Write current geometry
		GeoFile::AmberCrd.write(@geometry, :file => calc_dir+"/crd")

		# Forcedump.dat must be deleted, otherwise is not updated
		FileUtils.rm(calc_dir + "/forcedump.dat") if FileTest.exist?(calc_dir + "/forcedump.dat")

		# Run sander
		ld_path = ""
		ld_path = " export LD_LIBRARY_PATH=#{@settings[:amber_ld_path]}:$LD_LIBRARY_PATH;" if @settings[:amber_ld_path] != ""
		system("cd #{calc_dir}; export AMBERHOME=#{@settings[:amber_amberhome]};#{ld_path} #{@settings[:amber_sander_exe]} -O -i sander.in -p prmtop -o sander.out -c crd -r restrt > /dev/null 2> sander.err")
		Cuby::error("Calculation failed, sander did not produce the file sander.out",self) unless FileTest.exist?(calc_dir + "/sander.out")
		unless FileTest.exist?(calc_dir + "/forcedump.dat")
			amber_check_sander_out_errors
			Cuby::error("Calculation failed, sander did not produce the file forcedump.dat",self) 
		end
	end


	ENERGY_COMPONENTS = [
		:total_energy,
		:van_der_Waals,
		:electrostatics,
		:GBM,
		:bonds,
		:angles,
		:dihedrals,
		:"1-4_van_der_Waals",
		:"1-4_electrostatics"
	]

	def amber_read_results
		results = Results.new
		results.gradient = Gradient.new if @what.include?(:gradient)

		energyOK = false
		gradOK = false
		readgrad = false
		readenergy = false
		energycount = 0
		gradlinecount = 0
		IO.foreach(calc_dir + "/forcedump.dat") { |frcline|
			frcline.chomp!
			if frcline =~ /START of Energies/
				readenergy = true
				next
			end
			if readenergy
				if frcline =~ /[-]*[0-9]+[.][0-9]+E[+-][0-9]+/
					if energycount == 0
						results.energy = ($&).to_f
					else
						results.energy_components[ENERGY_COMPONENTS[energycount]] = ($&).to_f
						if energycount == 8
							energyOK = true
							readenergy = false
						end
					end
				end
				if  frcline =~ /NaN/
					@geometry.write_pdb(:file => calc_dir + "/failed_geo.pdb")
					Cuby::error "AMBER calculation failed, found NaN as a result (calculation #{@name})"
				end
				energycount += 1
			end
			if frcline =~ /Total Force/ && @what.include?(:gradient)
				readgrad = true
				next
			end
			if readgrad 
				gx, gy, gz = frcline.strip.split(/\s+/)  # split line, whitespace separated
				results.gradient[gradlinecount] = Coordinate.new(
					-1* gx.to_f,
					-1* gy.to_f,
					-1* gz.to_f
				)
				gradlinecount += 1
				if gradlinecount == @geometry.size
					gradOK = true
					break
				else
					gradOK = false
				end
			end
		}

		# Check if the results were read successfuly
		Cuby::error("Energy not found in file \"forcedump.dat\"", self) if ! energyOK
		Cuby::error("Valid gradient not found in file \"forcedump.dat\"", self) if ! gradOK && @what.include?(:gradient)

		# Calculate the interaction with point charges
		if @point_charges
			pch_energy = @point_charges.interaction_with_molecule(@geometry, @atomic_charges, results.gradient, @what.include?(:point_charges_gradient))
			results.energy_components[:point_charges_interaction] = pch_energy
			results.energy += pch_energy
		end

		return results
	end

	def amber_check_sander_out_errors
		parser = OutputParser.new(in_calc_dir("sander.out"), @name)
		parser.add_error(/Regenerate prmtop file with bondi radii/, "The calculation requires bondii radii in the prmtop file.")
		parser.execute
	end

	def amber_check_leap_out_errors
		parser = OutputParser.new(in_calc_dir("tleap.out"), @name)
		parser.add_error(/Could not open file/, "Error in parmtop creation, leap could not find some file\n(check file #{in_calc_dir('tleap.out')} for details)")
		parser.add_error(/Leap added 1 missing atom/, "Error in parmtop creation, leap added an atom\n(this means the geometry does not match the residue template in the ff.,\ncheck file #{in_calc_dir('tleap.out')} for details))")
		parser.add_error(/Leap added .* missing atoms/, "Error in parmtop creation, leap added multiple atoms\n(this means the geometry does not match the residue template in the ff.)")
		parser.execute
	end

	def amber_check_config
		@settings.check_dir(:amber_amberhome)
	end

	def amber_generate_prmtop
		# Write a PDB file 
		pdb_file_name = "tmp_in.pdb"
		@geometry.write_pdb(:file => in_calc_dir(pdb_file_name))

		# Check all leaprc files
		leaprc_files = []
		@settings[:amber_leaprc].each{|filename|
			filename = filename.gsub("$AMBERHOME", @settings[:amber_amberhome])
			filename.gsub!("%cuby",Cuby.install_dir)
			filename.gsub!("%interface",interface_dir)
			filename.gsub!("%amberhome",@settings[:amber_amberhome])
			filename.gsub!("%pwd",Dir.pwd)
			filename.gsub!("%home",File.expand_path("~"))
			filename = File.expand_path(filename)
			unless FileTest.exist?(filename)
				Cuby::error("Leaprc file not found:\n#{}")
			end
			leaprc_files << filename
		}

		# Create a tleap input file containing the commands needed to build the prmtop
		File.open(in_calc_dir('tleap.in'),'w+') {|f| 
			leaprc_files.each{|filename|
				f.puts "#==============================================================================="
				f.puts "# From #{filename}"
				File.open(filename, "r"){|orig|
					while line = orig.gets do
						line.chomp!
						line.gsub!("%cuby",Cuby.install_dir)
						line.gsub!("%interface",interface_dir)
						line.gsub!("%amberhome",@settings[:amber_amberhome])
						line.gsub!("$AMBERHOME",@settings[:amber_amberhome])
						line.gsub!("%pwd",Dir.pwd)
						line.gsub!("%home",File.expand_path("~"))
						f.puts line
					end
				}
				f.puts
			}
			f.puts "#==============================================================================="
			f.puts "# Prepare calculation"
			f.puts "x = loadpdb #{pdb_file_name}"
			#f.puts "set default OldPrmtopFormat on" if @settings [:amber_old_top]

			# IGB7 solvent model: set radii
			if @settings[:solvent_model] == :igb7
				f.puts "set default PBradii bondi"
			end

			f.puts "saveamberparm x prmtop restrt.tmp"
			f.puts "quit"
		}

		# Run tleap
		ld_path = ""
		ld_path = " export LD_LIBRARY_PATH=#{@settings[:amber_ld_path]}:$LD_LIBRARY_PATH;" if @settings[:amber_ld_path] != ""
		system("cd #{calc_dir}; export AMBERHOME=#{@settings[:amber_amberhome]};#{ld_path} #{@settings[:amber_amberhome]}/bin/tleap -f tleap.in > tleap.out 2> tleap.err")

		# Check if prmtop was crated
		if FileTest.exist?(in_calc_dir("prmtop"))
			if File.size(in_calc_dir("prmtop")) == 0
				amber_check_leap_out_errors
				Cuby::error("Amber prmtop is empty")
			end
		else
			Cuby::error("Amber prmtop file was not created.")
		end
		amber_check_leap_out_errors

		# Delete unnecessary files
		File.delete(in_calc_dir("restrt.tmp"))

		return nil
	end

	def amber_read_charges
		# Reads the atomic charges from the prmtop file
		charges = AtomicCharges.new(:forcefield, @geometry.size)

		top = Prmtop.new(in_calc_dir("prmtop"))
		Cuby::error "Can't load MM charges from prmtop - different number of atoms (system #{name})" if @geometry.size != top.natom

		@geometry.each_index { |i|
			charges[i] = top.charges[i]
		}

		return charges
	end

	#=======================================================================
	# Prmtop reader/writer
	#=======================================================================
	
	# It is now used only to get charges from the prmtop file, but can be used
	# to change them and some other properties

	class Prmtop
		#: Container for prmtop file. Reads some properties from the top, converting them
		#: to cuby units. These values can be modified and written back to the file.

		attr_reader :natom, :charges, :radii, :screens
		attr_writer :charges, :radii, :screens

		def initialize(filename)
			read_top(filename)
		end

		def read_top(topfilename)
			#: Read data from prmtop, filing @natom, @charges, @radii and @screens
			@toplines = IO.readlines(topfilename)

			# get pointers
			@pointers = []
			i = 0
			while @toplines[i] !~ /FLAG POINTERS/
				i += 1
			end
			i += 2
			while @toplines[i] !~ /%/
				items = @toplines[i].split
				items.each { |item| @pointers.push(item.to_i)}
				i += 1
			end
			@natom = @pointers[0]

			#get charges
			@charges = []
			while @toplines[i] !~ /FLAG CHARGE/
				i += 1
			end
			i += 2
			while @toplines[i] !~ /%/
				items = @toplines[i].split
				items.each { |item| @charges << ((item.to_f / 18.2223 * 100000000).round.to_f / 100000000)}
				i += 1
			end

			# get radii
			@radii = []
			while @toplines[i] !~ /FLAG RADII/
				i += 1
			end
			i += 2
			while @toplines[i] !~ /%/
				items = @toplines[i].split
				items.each { |item| @radii << item.to_f }
				i += 1
			end

			# get screens
			@screens = []
			while @toplines[i] !~ /FLAG SCREEN/
				i += 1
			end
			i += 2
			while (@toplines[i] !~ /%/) && (i < @toplines.length)
				items = @toplines[i].split
				items.each { |item| @screens << item.to_f }
				i += 1
			end
		end

		def write_top(topfilename)
			#: Write prmtop, overwriting the respective sections with data from
			#: @charges, @radii and @screens
			file = File.new(topfilename,"w+")
			i = 0
			while @toplines[i] !~ /FLAG CHARGE/
				file.puts @toplines[i]
				i += 1
			end
			file.puts @toplines[i]
			i += 1
			file.puts @toplines[i]
			i += 1
			#now on first line, write out new charges
			@charges.each_index { |chi|
				file.printf("%16.8E",@charges[chi] * 18.2223)
				if (chi+1)%5 == 0
					file.printf("\n")
				else
					file.printf("\n") if chi == @charges.length - 1
				end
			}
			#seek charges end
			while @toplines[i] !~ /%/
				i += 1
			end
			# seek radii start
			while @toplines[i] !~ /FLAG RADII/
				file.puts @toplines[i]
				i += 1
			end
			file.puts @toplines[i]
			i += 1
			file.puts @toplines[i]
			i += 1
			#now on first line, write out new radii
			@radii.each_index { |chi|
				file.printf("%16.8E",@radii[chi])
				if (chi+1)%5 == 0
					file.printf("\n")
				else
					file.printf("\n") if chi == @radii.length - 1
				end
			}
			#seek radii end
			while @toplines[i] !~ /%/
				i += 1
			end
			
			# seek screen start
			while @toplines[i] !~ /FLAG SCREEN/
				file.puts @toplines[i]
				i += 1
			end
			file.puts @toplines[i]
			i += 1
			file.puts @toplines[i]
			i += 1
			#now on first line, write out new screen
			@screens.each_index { |chi|
				file.printf("%16.8E",@screens[chi])
				if (chi+1)%5 == 0
					file.printf("\n")
				else
					file.printf("\n") if chi == @screens.length - 1
				end
			}
			#seek screen end
			while (@toplines[i] !~ /%/) && (i < @toplines.length)
				i += 1
			end
			
			#continue with old top
			while i < @toplines.length
				file.puts @toplines[i]
				i += 1
			end
			file.close
		end
	end
end
