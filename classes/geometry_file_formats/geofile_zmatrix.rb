################################################################################
#
# Module GeoFile::ZMatrixFile
#
# Author: Jan Rezac
# Date created: 2009-06-25
# License: Cuby license
# Description: ZMatrixFile file reader/writer
# Status: Working
#
################################################################################

require "classes/internal_coordinates/z_matrix"

module GeoFile
module ZMatrixFile
	def ZMatrixFile.read(geometry, arguments = {})
		#* :file - could be either File (file handle) or string (file name)
		#* :read_into_geometry
		#* :do_not_save_zmatrix - if set to true, parent zmatrix is not linked to the geometry

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

		s = file.read
		file.close if close
		zmat = ZMatrix.from_s(s)
		zmat.to_geometry.each{|atom| geometry << atom}

		if read_into
			geometry.coords_from_geometry(readgeo)
		end

		geometry.info[:z_matrix] = zmat unless arguments[:do_not_save_zmatrix]
	end
end
end
