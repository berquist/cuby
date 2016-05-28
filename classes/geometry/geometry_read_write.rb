################################################################################
#
# Module GeometryReadWrite
#
# Author: Jan Rezac
# Date created: 2008-11-03
# License: Cuby license
# Description: Geometry's routines for reading and writing files
# Extends: Geometry
# Status: Works, Documented
#
################################################################################

#: Module extending Geometry with methods to read and write the geometry in
#: common formats.

#: The file formats are separate modules within module GeoFile. They can be used
#: either directly, or via interface provided by GeoFile module. Here, these
#: methods are wrapped to provide direct access to file reading/writing from
#: Geometry class.

require "classes/geometry_file_formats/geofile.rb"

module GeometryReadWrite

	#=======================================================================
	# General read / write methods
	#=======================================================================
	
	def read_file(filename,arguments = {}) # => Symbol
		#: Reads geometry from file. Supports formats defined i GeoFile module, recognized by filename extension.
		#: Returns file type symbol.
		#
		#: Arguments are passed to the method reading the selected file format, see the
		#: details in corresponding file in file_formats directory.
		
		# Check if file exist
		if FileTest.exist?(filename)
			Cuby::error("File \"#{filename}\" is a directory") if FileTest.directory?(filename)
		else
			Cuby::error("File \"#{filename}\" does not exist")
		end

		# Get file type
		if arguments[:filetype]
			filetype = arguments[:filetype]
			arguments[:file] = filename
			GeoFile::read(self, filetype, arguments)
			@info[:filetype] = filetype
			return filetype
		else
			filetype = GeoFile.type_from_filename(filename)
			if filetype
				# Read the file from a file
				arguments[:file] = filename
				GeoFile::read(self, filetype, arguments)
				@info[:filetype] = filetype
				return filetype
			else
				filetype, pipe = GeoFile.type_autodetect(filename)
				Cuby::error "File format detection failed - can't read geometry. Hint: name it *.xyz or *.pdb" if filetype == nil
				# Read the file from a pipe provided by autodetect
				arguments[:file] = pipe
				GeoFile::read(self, filetype, arguments)
				pipe.close
				@info[:filetype] = filetype
				return filetype
			end
		end

	end

	def write_file(filename, filetype = :xyz, arguments = {})
		#: Wrapper for writting the geometry in a format secified by the filetype argument.
		#: Arguments are passed to the method writing the selected file format, see the
		#: details in corresponding file in file_formats directory.
		arguments[:file] = filename
		GeoFile::write(self, filetype, arguments)
	end

	#=======================================================================
	# Specific file readers/writers
	#=======================================================================
	
	#: Following methods use the file formats defined in GeoFile directly

	def write_xyz (arguments = {})
		return GeoFile::XYZ.write(self,arguments)
	end

	def read_xyz (arguments = {})
		return GeoFile::XYZ.read(self,arguments)
	end

	def read_pdb (arguments = {})
		return GeoFile::PDB.read(self, arguments)
	end

	def write_pdb (arguments = {})
		return GeoFile::PDB.write(self, arguments)
	end

	def read_turbomole_coord(arguments = {})
		return GeoFile::TurbomoleCoord.write(self, arguments)
	end

	def write_turbomole_coord(arguments = {})
		return GeoFile::TurbomoleCoord.write(self, arguments)
	end

	def write_amber_crd (arguments = {})
		return GeoFile::AmberCrd.write(self, arguments)
	end

end
