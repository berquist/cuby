################################################################################
#
# CHARMM interface
#
# Author: Jan Rezac
# Date created: 2014-01-14
# License: Cuby4 license
# Description: Interface for external calculations in CHARMM
# Status: Development
#
################################################################################

#===============================================================================
# Interface to the CHARMM molecular mechanics package
# Commercial software
# http://www.charmm.org/
#===============================================================================

require "erb" # ERB templating system, used for construction of the input
require "classes/tools/output_parser.rb"
require "classes/misc/search_for_file.rb"

# ToDo
# generate - steram files
# generate - extra input, separate for the generation
# input construction - gradient with more decimals
# test everything

module InterfaceCharmm
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :warning
	DEVELOPMENT_STATUS = "Interface works but the results were not tested"
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :gradient, :atomic_charges, :point_charges, :point_charges_gradient]
	MODIFIER = false
	DIRECTORY = "charmm"
	# Methods provided by the interface:
	METHODS = {
		:forcefield		=> [:energy, :gradient, :atomic_charges, :point_charges, :point_charges_gradient]
	}
	ATOMIC_CHARGES = [:forcefield]
	#=======================================================================

	CHARMM_INP_FN = "charmm.inp"
	CHARMM_OUT_FN = "charmm.out"
	
	def prepare_interface
		# Prepare calculation directory
		if calc_dir_mkdir(CHARMM_INP_FN, CHARMM_OUT_FN) == :old
			calc_using_input
			read_charmm_psf(calc_dir + "/sys.psf")
			return
		end

		charmm_write_input(what)

		if @settings.set?(:charmm_psf_file)
			# Copy PSF file
			if FileTest.exist?(@settings[:charmm_psf_file])
				FileUtils.cp(@settings[:charmm_psf_file], calc_dir + "/sys.psf")	
			else
				Cuby::error("CHARMM PSF file not found")
			end
		else
			# Genrate the PSF file automatically
			charmm_generate_psf
		end

		read_charmm_psf(calc_dir + "/sys.psf")

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		charmm_run_calculation
		
		results = charmm_read(what)

		if what.include?(:atomic_charges)
			results.atomic_charges = @atomic_charges
		end
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
		return @atomic_charges
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def charmm_stream_block(keyword)
		# Prepare a list of charmm stream files used
		stream_block = ""
		@settings[keyword].each{|fn|
			if FileTest.file?(File.expand_path(fn))
				fn = File.expand_path(fn)
			elsif fn = SearchForFile.search(:filename => fn, :locations => @settings[:charmm_ff_paths], :raise_error => false)
				fn = File.expand_path(fn)
			else
				Cuby::error("Charmm stream file #{fn} does not exist or is not a file")
			end
			stream_block << "stream \"#{fn}\"\n"
		}
		return stream_block
	end


	def top_par_files
		# Resolve names of top and par
		top = File.expand_path(@settings[:charmm_ff_top])
		par = File.expand_path(@settings[:charmm_ff_par])

		# Search for top
		unless FileTest.file?(top)
			top = SearchForFile.search(:filename => @settings[:charmm_ff_top], :locations => @settings[:charmm_ff_paths], :raise_error => false)
			Cuby::error("Charmm top file not found.") unless top
		end

		# Search for par
		unless FileTest.file?(par)
			par = SearchForFile.search(:filename => @settings[:charmm_ff_par], :locations => @settings[:charmm_ff_paths], :raise_error => false)
			Cuby::error("Charmm par file not found.") unless par
		end

		return [top, par]
	end

	def charmm_write_input(what)
		# Generate the input file from a template
		top, par = top_par_files
		template = ERB.new(IO.read(interface_dir + "/templates/charmm.inp.erb"))
		file = File.open(in_calc_dir(CHARMM_INP_FN),"w+")
		file.print(template.result(binding))
		file.close
	end
	

	def charmm_run_calculation
		# Write geometry
		write_charmm_coord(calc_dir + "/sys.crd")

		# Run calculation
		if @settings[:parallel] > 1
			# Parallel
			if ! system("cd #{calc_dir} && #{@settings[:charmm_exe_mpi]} #{@settings[:parallel]} #{CHARMM_INP_FN} 2> charmm.err > charmm.out")
				Cuby::error ("CHARMM returned nonzero exit code (#{@name})")
			end
		else
			# Serial
			unless system("cd #{calc_dir}; #{@settings[:charmm_exe]} < #{CHARMM_INP_FN} 2> charmm.err > #{CHARMM_OUT_FN}")
				Cuby::error "CHARMM returned nonzero exit code (#{@name})"
			end
		end

	end

	def charmm_read(what)
		# Initialize empty results
		results = Results.new
		results.gradient = Gradient.zero(@geometry.size) if @what.include?(:gradient)

		parser = OutputParser.new(in_calc_dir(CHARMM_OUT_FN), @name)

		parser.add_pattern(:energy, /^ENER> +0 *([-]*[0-9]+\.[0-9]+)/, :get => :submatch)

		if what.include?(:gradient)
			parser.add_block(:gradient_x, /scalar dx/, /^\s\(/, /CHARMM>\s*$/, {:count => @geometry.size})
			parser.add_block(:gradient_y, /scalar dy/, /^\s\(/, /CHARMM>\s*$/, {:count => @geometry.size})
			parser.add_block(:gradient_z, /scalar dz/, /^\s\(/, /CHARMM>\s*$/, {:count => @geometry.size})
		end

		parser.execute{|name, result, count|
			case name
			when :gradient_x
				results.gradient[count].x = result.split.last.to_f
			when :gradient_y
				results.gradient[count].y = result.split.last.to_f
			when :gradient_z
				results.gradient[count].z = result.split.last.to_f
			end
		}

		results.energy = parser[:energy]

		# Add interaction of with external point charges
		# Calculated by cuby using MM charges read from PSF
		if @point_charges
			pch_energy = @point_charges.interaction_with_molecule(@geometry, @atomic_charges, results.gradient, @what.include?(:point_charges_gradient))
			results.energy_components[:point_charges_interaction] = pch_energy
			results.energy += pch_energy
		end

		return results
	end

	def write_charmm_coord(filename)
		#example
		#    6    1 ALA  HA     5.64102  -7.88010  -4.79229 A    1      1.00000
		f = File.open(filename, "w+")
		# Write header
		f.puts "* NONE *"
		f.puts "*  DATE:    00/00/00     00:00:00      CREATED BY USER: cuby"
		f.puts "*"
		# Write # of atoms
		if @geometry.size < 100000
                	f.printf("%5d\n",@geometry.size)
			crd_line_format = "%5d%5d %-4s %-4s%10.5f%10.5f%10.5f %-4s %-4d%10.5f\n"
                else
			f.printf("%10d EXT\n",@geometry.size)
			crd_line_format = "%10d%10d  %-8s  %-8s%20.10f%20.10f%20.10f  %-8s  %-8d%20.10f\n"
                end 
		# Write coordinates
		recno = 0
		lastresno = 0
                lastsegid = "XXX"
                lastresname = "XXX"
		@geometry.each{|atom|
			recno += 1 if atom.properties[:pdb_res_no] != lastresno or atom.properties[:pdb_segid] != lastsegid or atom.properties[:pdb_res_name] != lastresname
			lastresno = atom.properties[:pdb_res_no]
                        lastsegid = atom.properties[:pdb_segid]
                        lastresname  = atom.properties[:pdb_res_name]
			f.printf(crd_line_format,
				atom.properties[:pdb_atom_no],
				recno,
				atom.properties[:pdb_res_name],
				atom.properties[:pdb_atom_name],
				atom.x,
				atom.y,
				atom.z,
				atom.properties[:pdb_segid],
				atom.properties[:pdb_res_no],
				1.0
			)
		}
		f.close
	end

	def read_charmm_psf(file)
		#: Read PSF file and save atom attributes, including its charge, in the geometry
		#: This information is later needed for building correct coordinate file for
		#: CHARMM.

		# Usual setup to use both File or filename as argument
		close = false
		f = file
		if f.class == String
			f = File.open(f,"r")
			close = true
		end

		# Atomic charges object
		@atomic_charges = AtomicCharges.new(:forcefield, @geometry.size)

		# Reading the file line by line
		atoms = false
		while line = f.gets
			line.chomp!

			# Skip everything after list of atoms
			break if atoms && line == ""

			# Read atom properties
			if atoms
				num, segid, resno, resname, atname, unknown, charge, mass, other = line.split
				num = num.to_i
				resno = resno.to_i
				charge = charge.to_f

				atom = @geometry[num - 1]
				atom.properties[:pdb_atom_no] = num
				atom.properties[:pdb_res_no] = resno
				atom.properties[:pdb_res_name] = resname
				atom.properties[:pdb_atom_name] = atname
				atom.properties[:pdb_segid] = segid
				@atomic_charges[num -1] = charge
			end

			# Search for start of atom list section
			if line =~ /NATOM/
				atoms = true
				natom = line.split[0].to_i
				if natom != @geometry.size
					Cuby::error "Number of atoms in PSF file does not agree with number of atoms in the geometry"
				end
			end
		end

		f.close if close
	end

	#=======================================================================
	# Private methods: Automated PSF generation
	#=======================================================================

	# Amino acids - when found on the ends of segment, the termination should be specified
	AMINO_ACIDS = [
		"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU",
		"GLY", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", 
		"MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR",
		"VAL"
	]

	def charmm_generate_psf
		# Check PDB data in geometry
		unless @geometry[0].properties[:pdb_atom_name]
			Cuby::error("Automated PSF generation requires a geometry in PDB format")
		end

		# Check segid, supply one if missing
		unless @geometry[0].properties[:pdb_segid]
			Cuby::warning("No segment ID found in the input PDB geometry, setting it to 'A' for all atoms")
			@geometry.each{|atom|
				atom.properties[:pdb_segid] = "A"
			}
		end

		# Write a PDB for the generation
		@geometry.write_pdb(:file => in_calc_dir('tmp_all.pdb'), :extra_columns => true)

		# Split the geometry into the PDB segments
		segment_names, segment_geo = split_pdb_to_segments

		# Write all segments
		segment_geo.each_index{|i|
			segment_geo[i].write_pdb(:file => in_calc_dir("tmp_#{i}.pdb"), :extra_columns => true)
		}

		# Generate the charmm input file from a template
		top, par = top_par_files
		template = ERB.new(IO.read(interface_dir + "/templates/generate_psf.inp.erb"))
		file = File.open(in_calc_dir('generate_psf.inp'),"w+")
		file.print(template.result(binding))
		file.close

		# Execute the generation in calculation directory
		unless system("cd #{calc_dir}; #{@settings[:charmm_exe]} < generate_psf.inp 2> generate_psf.err > generate_psf.out")
			Cuby::error "CHARMM returned nonzero exit code when generating PSF from PDB(#{@name})"
		end

		# Load the created PDB for check
		newpdb = Geometry.new
		newpdb.read_pdb(:file => in_calc_dir("generated.pdb"), :extra_columns => true)

		# Check number of atoms
		if newpdb.size > @geometry.size
			Cuby::error "CHARMM: Problem in generation PSF file,
#{newpdb.size - @geometry.size} more atoms generated than in original geometry"
		elsif newpdb.size < @geometry.size
			Cuby::error "CHARMM: Problem in generation PSF file,
#{-newpdb.size + @geometry.size} less atoms generated than in original geometry"
		end

		# Check matching elements
		@geometry.each_index{|i|
			unless @geometry[i].element == newpdb[i].element
				Cuby::error "CHARMM: Problem in generation PSF file, element does not match (atom #{i+1})"
			end
			unless @geometry[i].properties[:pdb_atom_name] == newpdb[i].properties[:pdb_atom_name]
				Cuby::error "CHARMM: Problem in generation PSF file, atom name does not match (atom #{i+1})"
			end
		}


	end

	def split_pdb_to_segments
		# Divide a geometry into separate geometries using the
		# PDB segment id information.
		# Returns array of segment names and array of geometries
		segment_names = []
		segment_geo = []
		i = 0
		segment_names << @geometry[0].properties[:pdb_segid]
		segment_geo[i] = Geometry.new
		segment_geo[i] << @geometry[0]
		(@geometry.size - 1).times{|i|
			if @geometry[i+1].properties[:pdb_segid] != segment_names.last
				segment_names << @geometry[i+1].properties[:pdb_segid]
				segment_geo << Geometry.new
			end
			segment_geo.last << @geometry[i+1]
		}
		Cuby::log.puts_debug("PDB segments found:")
		segment_names.each{|n| Cuby::log.puts_debug("   #{n}")}

		return [segment_names, segment_geo]
	end

	def segment_terminal_patches(geo, segname)
		# User input / autodetection of terminal patches for peptide chains
		# Used in the ERB template
		seg_start = "none"
		seg_end = "none"
		if @settings[:charmm_segment_end_patches].has_key?(segname)
			# User input
			Cuby::log.puts_debug("Terminal patches from input used for segment #{segname}")
			string = @settings[:charmm_segment_end_patches][segname]
			# Check format
			unless string =~ /^[a-z,A-Z,0-9]+\s*,\s*[a-z,A-Z,0-9]+$/
				Cuby::error "Invalid format of entry #{segname} in keyword charmm_segment_end_patches"
			end
			# Assign patches
			seg_start, seg_end = string.strip.split(/\s*,\s*/)
		else
			# Autodetection
			Cuby::log.puts_debug("Terminal patches autodetection for segment #{segname}")
			if AMINO_ACIDS.include?(geo[0].properties[:pdb_res_name])
				# Get first residue
				g = geo.geometry_from_selection(":1")
				if g.atomlist_from_selection("%atomname(HT1)").size == 1 &&
				   g.atomlist_from_selection("%atomname(HT2)").size == 1 &&
				   g.atomlist_from_selection("%atomname(HT3)").size == 1
					seg_start = 'NTER'
				elsif g.atomlist_from_selection("%atomname(HT1)").size == 1 &&
				   g.atomlist_from_selection("%atomname(HT2)").size == 1 &&
				   g.atomlist_from_selection("%atomname(HT3)").size == 0
					seg_start = 'NNEU'
				end
			end
			if AMINO_ACIDS.include?(geo.last.properties[:pdb_res_name])
				g = geo.geometry_from_selection(":#{geo.last.properties[:pdb_res_no]}")
				if g.atomlist_from_selection("%atomname(OT1)").size == 1 &&
				   g.atomlist_from_selection("%atomname(OT2)").size == 1 &&
				   g.atomlist_from_selection("%atomname(HT2)").size == 0
					seg_end = 'CTER'
				elsif g.atomlist_from_selection("%atomname(OT1)").size == 1 &&
				   g.atomlist_from_selection("%atomname(OT2)").size == 1 &&
				   g.atomlist_from_selection("%atomname(HT2)").size == 1
					seg_end = 'CNEU'
				end
			end
			Cuby::log.puts_debug(" using 'first #{seg_start} last #{seg_end}'")
		end
		return [seg_start, seg_end]
	end

	def segment_water_setup(geo, segname)
		if geo[0].properties[:pdb_res_name] == "TIP3"
			# TIP3P water
			# Check the segment is all water
			geo.each{|atom|
				unless atom.properties[:pdb_res_name] == "TIP3"
					Cuby::error "CHARMM PSF generation: water molecules must be in separate segment
(residue #{atom.properties[:pdb_res_name]} found in segment #{segname})"
				end
			}
			waterfix = 'noangle nodihedral'
		else
			# Non-water segment should not contain water
			waterfix = ''
			geo.each{|atom|
				if atom.properties[:pdb_res_name] == "TIP3"
					Cuby::error "CHARMM PSF generation: water molecules must be in separate segment
(residue #{atom.properties[:pdb_res_name]} found in segment #{segname})"
				end
			}
		end
		return waterfix
	end


end

