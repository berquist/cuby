################################################################################
#
# Module GeoFile::TurbomoleCoord
#
# Author: Jan Rezac
# Date created: 2009-06-25
# License: Cuby license
# Description: Turbomole coord file reader/writer
# Status: Working
#
################################################################################

module GeoFile
module TurbomoleCoord
	def TurbomoleCoord.read(geometry, arguments = {})
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

		# Check first line
		line = file.gets
		Cuby::error("File not in turbomole coord format") if line !~ /^\$coord/


		i = 0
		while line = file.gets do
			break if line =~ /\$/
			x,y,z,e = line.split
			readgeo << Atom.new(e,
					    x.to_f * BOHR2ANGSTROM, 
					    y.to_f * BOHR2ANGSTROM, 
					    z.to_f * BOHR2ANGSTROM,
					    {:file_index => i}
					   )
			i += 1
		end

		if read_into
			geometry.coords_from_geometry(readgeo)
		end
	end

	def TurbomoleCoord.write(geometry, arguments = {})
		arguments[:file] = $stdout if arguments[:file].nil?

		file = arguments[:file]

		close = false
		if file.class == String
			close = true
			mode="w+"
			file = File.new(file,mode)
		end
		unless file.kind_of?(IO)
			raise(TypeError,"File must be IO object")
		end
		file.puts "$coord"
		geometry.each {|atom|
			file.printf("%20.14f  %20.14f  %20.14f%7s\n",
				    atom.x * ANGSTROM2BOHR,
				    atom.y * ANGSTROM2BOHR,
				    atom.z * ANGSTROM2BOHR,
				    atom.element.to_s.downcase)
		}
		file.puts "$end"
		file.close if close
		return nil
	end
end
end
