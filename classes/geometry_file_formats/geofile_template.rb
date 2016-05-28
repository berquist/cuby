################################################################################
#
# Module GeoFile::Template
#
# Author: Jan Rezac
# Date created: 2009-06-25
# License: Cuby license
# Description: Template file reader/writer
# Status: Working
#
################################################################################

module GeoFile
module Template
	def Template.read(geometry, arguments = {})
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

		while line = file.gets do
			### Read the file here
		end

		if read_into
			geometry.coords_from_geometry(readgeo)
		end
	end

	def Template.write(geometry, arguments = {})
		#* :append (true/false), default is false
		#* :second_line - string, default is ""

		arguments[:file] = $stdout if arguments[:file].nil?
		arguments[:append] = false if arguments[:append].nil?

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
		geometry.each { |atom| 
			file.printf("\n", atom.element.to_s,atom.x,atom.y,atom.z) 
		}
		file.close if close
		return nil
	end
end
end
