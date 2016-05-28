################################################################################
#
# Module GeoFile::DftbGen
#
# Author: Jan Rezac
# Date created: 2009-06-25
# License: Cuby license
# Description: Turbomole coord file reader/writer
# Status: Working
#
################################################################################

module GeoFile
module DftbGen

	def DftbGen.write(geometry, arguments = {})
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

		# Output starts here

		# get elements in structure - array of atoms
		presentatoms = geometry.elements_in_system

		# write file
		file.puts "#{geometry.size}  C"
		presentatoms.each { |a|
			file.printf("%s ",a.to_s.upcase)
		}
		file.puts
		geometry.each_index { |i|
			file.printf("%d %d %f %f %f\n",i+1,presentatoms.index(geometry[i].element)+1,geometry[i].x,geometry[i].y,geometry[i].z)
		}

		# Output ends here

		file.close if close
		return nil
	end
end
end
