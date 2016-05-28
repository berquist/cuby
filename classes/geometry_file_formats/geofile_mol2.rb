################################################################################
#
# Module GeoFile::MOL2
#
# Author: Jan Rezac
# Date created: 2009-08-17
# License: Cuby license
# Description: Tripos .mol2 file reader
# Status: Working
#
################################################################################

#: The .mol2 reader is limited ro reading the xyz coordinates from the file,
#: it ignores the connectivity and other information.

#: Format documentation: http://tripos.com/index.php?family=modules,SimplePage,,,&page=sup_mol2&s=0

module GeoFile
module MOL2
	def MOL2.read(geometry, arguments = {})
		#: Reads .mol2 file. Arguments:
		#* :file - could be either File (file handle) or string (file name)
		#* :read_into_geometry

		file = arguments[:file]
		read_into = arguments[:read_into_geometry]

		close = false
		if file.class == String
			close = true
			file = File.new(file,"r")
		end

		unless file.kind_of?(IO) || file.kind_of?(StringIO)
			raise(TypeError,"File must be IO object (is #{file.class})")
		end

		if read_into
			readgeo = Geometry.new
		else
			readgeo = geometry
			geometry.clear # delete old atoms
		end

		header_natom = 0
		# Seek the MOLECULE section
		while line = file.gets
			line.chomp!
			if line == "@<TRIPOS>MOLECULE"
				# Read name
				file.gets
				# Read number of atoms
				header_natom = file.gets.split[0].to_i
				break
			end
		end
		Cuby::error "Invalid MOL2 file - section MOLECULE not found" if line == nil
		Cuby::error "Invalid MOL2 file - no atoms" if header_natom == 0

		# Seek ATOM section
		while line = file.gets
			line.chomp!
			if line == "@<TRIPOS>ATOM"
				break
			end
		end
		Cuby::error "Invalid MOL2 file - section ATOM not found" if line == nil

		# Read the geometry
		total_charge = 0.0
		header_natom.times{|i|
			if (line = file.gets) == nil
				Cuby::error('Invalid MOL2 file - it contains less atoms than specified in the header')
			end
			line.chomp!
			id, name, x, y, z, type, sub_no, sub_name, charge, status = line.chomp.strip.split
			#element = type.sub(/\..*/,'')
			element = name.sub(/[0-9]+.*/,'')
			atom = Atom.new(element, x.to_f, y.to_f, z.to_f, {:file_index => i})
			atom.properties[:pdb_atom_name] = name
			atom.properties[:pdb_res_name] = sub_name
			atom.properties[:pdb_res_no] = sub_no
			charge = charge.to_f
			atom.properties[:charge] = charge
			readgeo << atom

			# Count charge
			total_charge += charge
		}
		geometry.info[:geofile_charge] = total_charge

		if read_into
			geometry.coords_from_geometry(readgeo)
		end
		return true
	end

	def MOL2.write(geometry, arguments = {})
		#: Writes geometry in MOL2 format. Arguments:
		#* :file - file handle or name, default is $stdout
		#* :bondlist - list of bonds
		##* :append

		arguments[:file] = $stdout if arguments[:file].nil?
		long = false
		long = arguments[:long] unless arguments[:long].nil?

		file = arguments[:file]

		# Deafults
		arguments[:name] = "UNK" unless arguments[:name]

		# Make bonds list
		raise "Bond list must be provided to mol2 writer" unless  arguments[:bondlist]
		bonds = arguments[:bondlist]

		# Open the file
		close = false
		if file.class == String
			close = true
			mode="w+"
			mode = "a+" if arguments[:append] == true
			file = File.new(file,mode)
		end
		unless file.kind_of?(IO)
			raise(TypeError,"File must be IO object")
		end

		# Molecule record
		file.puts "@<TRIPOS>MOLECULE"
		file.puts arguments[:name]
		file.printf("%5d %5d %5d %5d %5d\n",geometry.size, bonds.size, 1,0,0)
		file.puts "SMALL"
		file.puts "USER_CHARGES"
		file.puts
		file.puts

		# Name atoms
		#unless geometry[0].properties[:pdb_atom_name]
		#	geometry.pdb_name_atoms_by_num
		#end
		geometry.pdb_name_atoms_by_num

		# Atoms record
		file.puts "@<TRIPOS>ATOM"
		geometry.each_index{|i|
			atom = geometry[i]
			file.printf("%7d ",i+1) # Atom index
			file.printf("%-8s", atom.properties[:pdb_atom_name] ) # Atom name
			file.printf("%10.4f",atom.x)
			file.printf("%10.4f",atom.y)
			file.printf("%10.4f",atom.z)
			file.printf(" %-4s ",atom.element.to_s) # Atom type
			file.printf("%6d ",1)	# Substructure no.
			file.printf("%-6s","UNK") # Substructure name
			if atom.properties[:charge] == nil
				file.printf("%10.4f",0.0) # Charge
			else
				file.printf("%10.4f",atom.properties[:charge]) # Charge
			end
			file.printf("\n")
		}

		# Bonds
		file.puts "@<TRIPOS>BOND"
		bonds.each_index{|i|
			bond = bonds[i]
			file.printf("%6d",i+1)
			file.printf("%5d",bond[0]+1)
			file.printf("%5d",bond[1]+1)
			file.printf(" %-6s","1")
			file.printf("\n")
		}

		# Substructure
		file.puts "@<TRIPOS>SUBSTRUCTURE"
		file.puts "     1 UNK         1 TEMP              0 ****  ****    0 ROOT"



		file.close if close
		return nil
	end

end
end
