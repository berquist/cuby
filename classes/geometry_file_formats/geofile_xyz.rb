################################################################################
#
# Module GeoFile::XYZ
#
# Author: Jan Rezac
# Date created: 2009-06-25
# License: Cuby license
# Description: XYZ file reader/writer
# Status: Working
#
################################################################################

module GeoFile
module XYZ
	def XYZ.read(geometry, arguments = {})
		#: Reads XYZ file. Arguments:
		#* :file - could be either File (file handle) or string (file name)
		#* :read_into_geometry
		#* :velocities

		file = arguments[:file]
		read_into = arguments[:read_into_geometry]
		velocities = arguments[:velocities]

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

		# Try to read the file
		firstline = file.gets
		return false if firstline == nil

		# Try to get number of atoms
		firstline = firstline.chomp.strip
		if firstline =~ /^[0-9]+$/
			header_size = firstline.to_i
		elsif firstline =~ /\s*[a-zA-Z]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+/
			# No header!
			Cuby::warning("The file seems to be an .xyz file without a header. Reading until end of file or a non-atom line.")
			a = firstline.split
			atom = Atom.new(a[0],a[1].to_f,a[2].to_f,a[3].to_f)
			readgeo << atom
			while (line = file.gets) != nil
				if line =~ /\s*[a-zA-Z]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+/
					a = line.split
					atom = Atom.new(a[0],a[1].to_f,a[2].to_f,a[3].to_f)
					readgeo << atom
				else
					break
				end
			end
			if read_into
				geometry.coords_from_geometry(readgeo)
			end
			return true
		else
			Cuby::error('Invalid XYZ file - first line is not a number')
		end

		# Try to read second line
		if (secondline = file.gets) == nil
			Cuby::error('Invalid XYZ file - has only one line')
		end

		if secondline =~ /^\s*[a-z,A-Z,0-9]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+/
			Cuby::warning("Second line of the .xyz file looks like coordinates of an atom")
		end

		# Read the geometry
		header_size.times{|i|
			if (line = file.gets) == nil
				Cuby::error("Invalid XYZ file - it contains less atoms than specified in the header (#{header_size})")
			end
			line.chomp!
			if line =~ /\s*[a-zA-Z]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+/
				a = line.split
				atom = Atom.new(a[0],a[1].to_f,a[2].to_f,a[3].to_f)
				if velocities
					Cuby::error "Velocities expected in XYZ file, but only #{a.size - 4} extra columns found" unless a.size == 7
					atom.properties[:velocity] = Coordinate.new(a[4].to_f, a[5].to_f, a[6].to_f)
				end
				readgeo << atom
			elsif line =~ /\s*[a-zA-Z]+[0-9]*\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+/
				a = line.split
				e = a[0].gsub(/[0-9]/,"")
				atom = Atom.new(e,a[1].to_f,a[2].to_f,a[3].to_f)
				if velocities
					Cuby::error "Velocities expected in XYZ file, but only #{a.size - 4} extra columns found" unless a.size == 7
					atom.properties[:velocity] = Coordinate.new(a[4].to_f, a[5].to_f, a[6].to_f)
				end
				readgeo << atom
			elsif line =~ /\s*[0-9]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+\s+[-]?[0-9]+\.[0-9]+/
				# Number instead of element
				a = line.split
				e = PeriodicTable::ELEMENTS[a[0].to_i]
				atom = Atom.new(e,a[1].to_f,a[2].to_f,a[3].to_f)
				if velocities
					Cuby::error "Velocities expected in XYZ file, but only #{a.size - 4} extra columns found" unless a.size == 7
					atom.properties[:velocity] = Coordinate.new(a[4].to_f, a[5].to_f, a[6].to_f)
				end
				readgeo << atom
			else
				Cuby::error("Invalid XYZ file - wrong format of line:\n'#{line}'")
			end

		}

		# file_index numbering
		geometry.each_index{|i| geometry[i].properties[:file_index] = i}

		if read_into
			geometry.coords_from_geometry(readgeo)
			if velocities
				geometry.each_index{|i| geometry[i].properties[:velocity] = readgeo[i].properties[:velocity]}
			end
		end

		geometry.info[:xyz_remark] = secondline.chomp.strip

		return true
	end

	def XYZ.write(geometry, arguments = {})
		#: Writes geometry in xyz format. Arguments:
		#* :file - file handle or name, default is $stdout
		#* :append (true/false), default is false
		#* :second_line - string, default is ""
		#* :format (string) - printf format used to print coordinate, default " %13.9f"
		#* :no_header (true/false), use to supress printing first two lines, default is false
		#* :velocities - write velocities (array of velocities as coordinates is passed as an argument)

		arguments[:file] = $stdout if arguments[:file].nil?
		arguments[:append] = false if arguments[:append].nil?
		arguments[:second_line] = "" if arguments[:second_line].nil?
		arguments[:format] = " %13.9f" if arguments[:format].nil?
		velocities = arguments[:velocities]

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
		unless arguments[:no_header] == true
			file.puts geometry.size
			file.puts arguments[:second_line]
		end
		geometry.each_with_index { |atom, i| 
			file.printf("%3s #{arguments[:format]}#{arguments[:format]}#{arguments[:format]}",
						   atom.element.to_s,atom.x,atom.y,atom.z) 
			if velocities
				file.printf("#{arguments[:format]}#{arguments[:format]}#{arguments[:format]}",
					    velocities[i].x,
					    velocities[i].y,
					    velocities[i].z
					   )
			end
			file.puts
		}
		file.close if close
		return nil
	end
end
end
