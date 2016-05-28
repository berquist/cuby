################################################################################
#
# Module GeoFile::AmberCrd
#
# Author: Jan Rezac
# Date created: 2009-06-25
# License: Cuby license
# Description: AMBER .crd file writer
# Status: Working
#
################################################################################

module GeoFile
module AmberCrd
	def AmberCrd.write(geometry, arguments = {})
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

		file.printf("\n%6d\n",geometry.size)
		eo = false
		geometry.each { |atom| 
			file.printf("%12.7f%12.7f%12.7f",atom.x,atom.y,atom.z)
			file.printf("\n") if eo 
			eo = !eo
		}
		file.printf("\n") if eo 
		file.close if close
		return nil
	end
end
end
