################################################################################
#
# Module GeoFile::SDF
#
# Author: Jan Rezac
# Date created: 2009-08-17
# License: Cuby license
# Description: SDF file reader
# Status: Working
#
################################################################################

#: The SDF reader is limited ro reading the xyz coordinates from the file,
#: it ignores the connectivity and other information.

module GeoFile
module SDF
	def SDF.read(geometry, arguments = {})
		#: Reads SDF file. Arguments:
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

		# Read the 3 comments
		3.times {
			firstline = file.gets
			return false if firstline == nil
		}

		# Read the information line
		firstline = file.gets
		return false if firstline == nil
		info = firstline.chomp.strip.split

		# Try to get number of atoms
		if info[0] =~ /^[0-9]+$/
			if info[0].size > 3
				header_natom = info[0][0..2].to_i
			else
				header_natom = info[0].to_i
			end
		else
			Cuby::error('Invalid SDF file - can''t read number of atoms record')
		end

		# Read the geometry
		header_natom.times{|i|
			if (line = file.gets) == nil
				Cuby::error('Invalid SDF file - it contains less atoms than specified in the header')
			end
			line.chomp!
			if matchdata =  / *([-]?[0-9]+\.[0-9]+) +([-]?[0-9]+\.[0-9]+) +([-]?[0-9]+\.[0-9]+) +([a-zA-Z]+) +.*/.match(line)
				readgeo << Atom.new(matchdata[4],matchdata[1].to_f, matchdata[2].to_f, matchdata[3].to_f, {:file_index => i})
			else
				Cuby::error("Invalid SDF file - wrong format of line:\n'#{line}'")
			end

		}

		if read_into
			geometry.coords_from_geometry(readgeo)
		end
		return true
	end

end
end
