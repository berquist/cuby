################################################################################
#
# Module GeoFile::PDB
#
# Author: Jan Rezac
# Date created: 2009-06-25
# License: Cuby license
# Description: PDB file reader/writer
# Status: Working
#
################################################################################

# This module provides read/write methods for PDB files

# http://www.wwpdb.org/documentation/format32/sect9.html
# http://bmerc-www.bu.edu/needle-doc/latest/atom-format.html

module GeoFile

module LPDB
	#: This module handles writting Long-coordinate PDBs, using filetype :lpdb
	
	def LPDB.read(geometry, arguments = {})
		return GeoFile::PDB.read(geometry, arguments)
	end

	def LPDB.write(geometry, arguments = {})
		arguments_copy = arguments.clone
		arguments_copy[:long] = true
		return GeoFile::PDB.write(geometry, arguments_copy)
	end
end

module PDB
	#LPDB_FORMAT = "%11.6f"
	#LPDB_LENGTH = 11
	
	# Write 15 decimals, read both
	LPDB_FORMAT = "%20.15f"
	LPDB_LENGTH = 20
	LPDB_LENGTH_OLD = 11

	def PDB.read(geometry, arguments = {})
		#: Long format (coordinates in %11.6f format) is autodetected
		#: arguments:
		#* :file (handle or filename)
		#* :element_heuristic # :auto, true, false
		#* :extra_columns (read extra column with additional info)
		#* :read_into_geometry
		#* :altloc (pick only one alternative location: nil - use A and print warning, letter - use selected, "*" - save both)

		heuristic = :auto
		heuristic = arguments[:element_heuristic] if arguments.has_key?(:element_heuristic)
		if arguments.has_key?(:extra_columns)
			extra_cols = arguments[:extra_columns]
		else
			extra_cols = true
		end
		read_into = arguments[:read_into_geometry]

		altloc = nil
		altloc = arguments[:altloc] if arguments.has_key?(:altloc)

		# create file from filename if necessary
		pdbfile = arguments[:file]
		close = false
		if pdbfile.class == String
			close = true
			begin	
				pdbfile = File.new(pdbfile,"r")
			rescue
				Cuby::error("Can't open PDB file \"#{pdbfile}\"")
				exit
			end
		end
		unless pdbfile.kind_of?(IO) || pdbfile.kind_of?(StringIO)
			raise(TypeError,"File must be IO object")
		end

		if read_into
			readgeo = Geometry.new
		else
			readgeo = geometry
			geometry.clear # delete old atoms
		end

		counter = 0
		alt_loc_count = 0
		linecount = 1
		uses_model = false
		end_found = false
		while line = pdbfile.gets do
			line.chomp!
			next if line == ''
			unless matchdata = /^(\S+)\s*/.match(line)
				Cuby::error("PDB format error at line #{linecount}")
			end
			linetype = matchdata[1]
			atom = nil
			case linetype
			when 'HETATM'
				Cuby::error("PDB format error at line #{linecount}: data after END") if end_found
				atom = PDB.line_to_atom(line, linecount, extra_cols, heuristic)
			when 'ATOM'
				Cuby::error("PDB format error at line #{linecount}: data after END") if end_found
				atom = PDB.line_to_atom(line, linecount, extra_cols, heuristic)
			when 'TER'
				if counter > 0
					geometry[counter-1].properties[:terafter] = true 
					#--------------------------------------------------
					# matchdata = /(.{6})(.{5})      (.{3}) (.)(.{4})./.match(line.ljust(27))
					# ternum = matchdata[2].to_i
					# resname = matchdata[3].strip
					# chain = matchdata[4]
					# resno = matchdata[5].to_i
					#-------------------------------------------------- 
				end
			when "MODEL"
				uses_model = true;
			when "REMARK"
				# header can be read here...
				# Read charge and multiplicity data
				rem, name, value = line.strip.split
				case name.downcase
				when "charge"
					geometry.info[:geofile_charge] = value.to_i
				when "multiplicity"
					geometry.info[:geofile_multiplicity] = value.to_i
				end
			when "END"
				if uses_model
					end_found = true
				else
					break 
				end
			when "ENDMDL"
				break
			end

			# Add atom if the line is ATOM or HETATM
			if atom
				# Handle alternate locations
				add_atom = false
				if atom.properties[:pdb_alt_loc]
					case altloc
					when nil
						if atom.properties[:pdb_alt_loc] == "A"
							add_atom = true 
							atom.properties[:pdb_alt_loc] = nil
							alt_loc_count += 1
						end
					when "A"
						add_atom = true if atom.properties[:pdb_alt_loc] == "A"
					when "B"
						add_atom = true if atom.properties[:pdb_alt_loc] == "B"
					when "*"
						add_atom = true
					end
				else
					add_atom = true
				end

				if add_atom
					readgeo[counter] = atom
					counter += 1
				end
			end
			linecount += 1
		end

		if alt_loc_count > 0
			Cuby::warning("PDB contained alternate location for #{alt_loc_count} atoms, variant A was picked by default")
		end

		# cloese file if opened here
		pdbfile.close if close

		# file_index numbering
		geometry.each_index{|i| geometry[i].properties[:file_index] = i}

		if read_into
			geometry.coords_from_geometry(readgeo)
		end

		# check number of atoms
		return false if counter == 0
		return true
	end

	def PDB.write(geometry, arguments = {})
		#: Writes geometry in PDB format. Arguments:
		#* :file - file handle or name, default is $stdout
		#* :extra_columns
		#* :long - use long coordinates
		#* :append
		#* :model_no
		#* :connectivity

		arguments[:file] = $stdout if arguments[:file].nil?
		long = false
		long = arguments[:long] unless arguments[:long].nil?

		file = arguments[:file]

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

		write_extended = arguments[:extra_columns]

		file.printf("MODEL%9d\n", arguments[:model_no]) if arguments[:model_no]

		counter = 1
		geometry.each { |atom| 
			file.puts PDB.atom_to_s(atom, write_extended, counter, long)
			if (atom.properties[:terafter]) #if atom is at end of residue
				file.puts("TER") 
			end
			counter += 1
		}

		if arguments[:connectivity]
			conn = arguments[:connectivity]
			geometry.each_index { |i| 
				atom = geometry[i]
				conn.bound_atoms(i).each_index{|j|
					if j % 4 == 0
						if j != 0
							file.puts
						end
						if atom.properties[:pdb_atom_no]
							atom_no = atom.properties[:pdb_atom_no]
						else
							atom_no = i + 1
						end
						file.printf("CONECT%5d", atom_no)
					end
					# Write bound atom numbers
					atom2 = geometry.at(conn.bound_atoms(i)[j])
					if atom2.properties[:pdb_atom_no]
						atom_no2 = atom2.properties[:pdb_atom_no]
					else
						atom_no2 = conn.bound_atoms(i)[j] + 1
					end
					file.printf("%5d", atom_no2)
				}
				if conn.bound_atoms(i).size > 0
					file.puts
				end
			}
		end

		if arguments[:model_no]
			file.puts("ENDMDL")
		else
			file.puts("END")
		end
		file.close if close
		return nil
	end

	PDB_ELEMENTS = {
		"x" => :X,
		"fe" => :Fe,
		"sod" => :Na,
		"na+" => :Na,
		"ip" => :Na,
		"f" => :F,
		"mg" => :Mg,
		"cl" => :Cl,
		"im" => :Cl,
		"zn" => :Zn,
		"mn" => :Mn,
		"br" => :Br
	}
	
	def PDB.element_from_type(string, line_no) # => element
		#: An heuristic procedure for determining element from PDB atom type.
		#: It uses a set of empiricaly derived rules and can fail when it
		#: encounters some unexpected atom type.
		
		element = string.downcase

		if element =~ /^[a-z]+[0-9]+[a-z]+/
			element.gsub!(/[0-9]+[a-z]+/,'')
		end

		element.gsub!(/[0-9,-]/,"") # remove numbers

		# Known elements lookup
		if PDB_ELEMENTS.has_key?(element)
			return PDB_ELEMENTS[element]
		end

		# remove non-alphabetic characters and some letters
		element.gsub!(/[*,',+,t,l,a,d,e,g,z,x,m,y]/,"")
		element.gsub!(/(.)w/,'\1') #remove all w's except first

		# remove label characters for HCNO
		first_char = element[0..0]
		if first_char == "o"
			element.gsub!(/o[h,p,q,c,r]+/,"o")
		elsif first_char == "h"
			element.gsub!(/h[h,n,c,b,o,p,q,r,y]+/,"h")
		elsif first_char == "n"
			element.gsub!(/n[h,c]+/,"n")
		elsif first_char == "c"
			element.gsub!(/c[b,h,q,c]+/,"c")
		end
		
		e =  element.downcase.capitalize.to_sym
		begin
			Atom.check_element(e)
		rescue
			Cuby::error("Can't get element from PDB atom type \"#{string}\" (line #{line_no})")
		end
		return e
	end

	def PDB.line_to_atom(string, line_no, extra_cols, element_heuristic)

		# Long format autodetection
		long = string[35..40] =~ /^[0-9]+$/
		xlong = string[35..49] =~ /^[0-9]+$/
		if xlong
			string = string.chomp.ljust(116)
			matchdata = /^(.{6})(.{5}) (.{4})(.{1})(.{4})(.{1})(.{5}).  (.{#{LPDB_LENGTH}})(.{#{LPDB_LENGTH}})(.{#{LPDB_LENGTH}})(.{6})(.{6}).{6}(.{4})(.{2})(.{2})/.match(string)
		elsif long
			string = string.chomp.ljust(89)
			matchdata = /^(.{6})(.{5}) (.{4})(.{1})(.{4})(.{1})(.{5}).  (.{#{LPDB_LENGTH_OLD}})(.{#{LPDB_LENGTH_OLD}})(.{#{LPDB_LENGTH_OLD}})(.{6})(.{6}).{6}(.{4})(.{2})(.{2})/.match(string)
		elsif
			# Standard format
			string = string.chomp.ljust(80)
			matchdata = /^(.{6})(.{5}) (.{4})(.{1})(.{4})(.{1})(.{5}).  (.{8})(.{8})(.{8})(.{6})(.{6}).{6}(.{4})(.{2})(.{2})/.match(string)
			if matchdata == nil
				# Try PDB with more than 1,000,000 atoms
				matchdata = /^(.{6})(.{6})(.{4})(.{1})(.{4})(.{1})(.{5}).  (.{8})(.{8})(.{8})(.{6})(.{6}).{6}(.{4})(.{2})(.{2})/.match(string)
			end
		end

		if matchdata == nil
			Cuby::error("PDB line #{line_no} can not be parsed:\n#{string}")
		end
	
		if matchdata[8] !~ /[0-9]/ || matchdata[9] !~ /[0-9]/ || matchdata[10] !~ /[0-9]/
			Cuby::error("Missing coordinates at PDB line #{line_no}")
		end
		x = matchdata[8].to_f
		y = matchdata[9].to_f
		z = matchdata[10].to_f

		case element_heuristic
		when true
			heuristic = true
		when false
			heuristic = false
		when :auto
			if matchdata[14] !~ /^\s*$/
				heuristic = false
			else
				 heuristic = true
			end
		else
			Cuby::error('Wrong value of parameter element_heuristic')
		end

		# Get element
		if heuristic
			element = PDB.element_from_type(matchdata[3].strip, line_no)
		else
			# Use element column or rigorous interpretation of atom name
			if matchdata[14] !~ /^\s*$/
				element = Atom.element_from_string(matchdata[14].strip)
			else
				es = /^(..)/.match(matchdata[3])[1]
				es.gsub!(/[0-9, ]/,'')
				begin
					element = Atom.element_from_string(es)
				rescue
					# Exceptions can be handled here
					element = PDB.element_from_type(matchdata[3].strip, line_no)
				end
			end
		end

		# Build properties
		properties = {
			:pdb_atom_no => matchdata[2].to_i,
			:pdb_atom_name => matchdata[3].strip,
			:pdb_res_name => matchdata[5].strip,
			:pdb_res_no => matchdata[7].to_i,
		}

		properties[:pdb_heteroatom] = true if matchdata[1] == 'HETATM'
		properties[:pdb_chain] = matchdata[6] if matchdata[6] !~ /^\s*$/
		properties[:pdb_alt_loc] = matchdata[4] if matchdata[4] !~ /\s/

		# Optional records
		if extra_cols
			properties[:pdb_occupancy] = matchdata[11].to_f if matchdata[11] !~ /^\s*$/
			properties[:pdb_temp_factor] = matchdata[12].to_f if matchdata[12] !~ /^\s*$/
			properties[:pdb_segid] = matchdata[13] if matchdata[13] !~ /^\s*$/
			properties[:pdb_element] = matchdata[14].strip if matchdata[14] !~ /^\s*$/
			properties[:pdb_charge] = matchdata[15] if matchdata[15] !~ /^\s*$/
		end

		return Atom.new(element, x,y,z, properties)
	end

	def PDB.atom_to_s(atom, write_extended, atom_no, long)
		# A modfication that allows residue names with four characters.
		# This does not conform to the PDB standard but the column occupied
		# by the fourth character is defined as empty, what is true when
		# proper residue names are used. Therefore, this option is always on.
		long_resname = true

		# Record type
		rectype = "ATOM"
		rectype = "HETATM" if atom.properties[:pdb_heteroatom]

		# Atom number
		if atom.properties[:pdb_atom_no]
			atom_no = atom.properties[:pdb_atom_no] 
		else
			raise "Atom number must be specified as argument if 'properties[:pdb_atom_no]' is mising" if atom_no == nil
		end

		# Atom name
		if atom.properties[:pdb_atom_name]
			if atom.element.to_s.size == 1 && atom.properties[:pdb_atom_name].size < 4 && atom.properties[:pdb_atom_name][0..0] !~ /[0-9]/
				atomname = " " + atom.properties[:pdb_atom_name]
			else
				atomname = atom.properties[:pdb_atom_name]
			end
		else
			if atom.element.to_s.size == 1
				atomname = " " + atom.element.to_s
			else
				atomname = atom.element.to_s
			end
		end

		# Residue name
		if atom.properties[:pdb_res_name]
			res_name = atom.properties[:pdb_res_name]
		else
			res_name = 'UNK'
		end

		if long_resname
			res_name = res_name[0..4] if res_name.size > 4
		else
			res_name = res_name[0..2] if res_name.size > 3
		end

		# Residue number
		if atom.properties[:pdb_res_no]
			res_no = atom.properties[:pdb_res_no]
		else
			res_no = 1
		end

		# Build format string
		# Atom number length
		if atom_no <= 99999
			# Normal
			format = "%-6s%5d "
		else
			# Extra decimal place (non-standard)
			format = "%-6s%6d"
		end
		# Residue name length
		if long_resname
			format += "%-4s%1s%-4s%1s"
		else
			format += "%-4s%1s%-3s %1s"
		end
		# Residue number
		if res_no <= 9999
			format += "%4d "
		else
			format += "%5d"
		end
		# Coordinates
		if long
			format += "   #{LPDB_FORMAT * 3}"
		else
			format += "   %8.3f%8.3f%8.3f"
		end

		# Write string
		s = sprintf(format,
			rectype,
			atom_no,
			atomname,
			atom.properties[:pdb_alt_loc],
			res_name,
			atom.properties[:pdb_chain],
			res_no,
			atom.x,
			atom.y,
			atom.z
		)
		# Extended info
		if write_extended
			if atom.properties[:pdb_occupancy]
				occupancy = atom.properties[:pdb_occupancy]
			else
				occupancy = 1.0
			end

			if atom.properties[:pdb_temp_factor]
				temp_factor = atom.properties[:pdb_temp_factor]
			else
				temp_factor = 0.0
			end

			s << sprintf("%6.2f%6.2f      %-4s%2s%2s",
				occupancy,
				temp_factor,
				atom.properties[:pdb_segid],
				atom.element.to_s,
				atom.properties[:pdb_charge]
			)
		end
		return s
	end

end
end
