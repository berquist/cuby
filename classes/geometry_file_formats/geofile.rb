################################################################################
#
# Module GeoFile
#
# Author: Jan Rezac
# Date created: 2009-06-25
# License: Cuby license
# Description: Namespace for molecular geometry file formats
# Status: Working, Documented
#
################################################################################


#: File reading methods return:
#* True when everything is OK
#* False when nothing was read
#: or raise an exception when reading error was encountered

require "stringio"

module GeoFile

	#: File types that can be read - filetype (symbol) => module name (string)
	TYPE_MODULE_R = {
		:pdb => "PDB",
		:lpdb => "LPDB",
		:xyz => "XYZ",
		:tm_coord => "TurbomoleCoord",
		:zmat => "ZMatrixFile",
		:sdf => "SDF",
		:mol2 => "MOL2"
	}

	#: File types that can be written - filetype (symbol) => module name (string)
	TYPE_MODULE_W = {
		:pdb => "PDB",
		:lpdb => "LPDB",
		:xyz => "XYZ",
		:tm_coord => "TurbomoleCoord",
		:amber_crd => "AmberCrd",
		:mol2 => "MOL2",
		:dftb_gen => "DftbGen"
	}

	#: Pattrens used to evaluate file type (symbol) from filename (string)
	TYPE_FN_PATTERNS = {
		:pdb => /.*\.pdb$/,
		:lpdb => /.*\.lpdb$/,
		:xyz => /.*\.xyz$/,
		:tm_coord => /^coord$|.*\.coord$/,
		:amber_crd => /^crd$|^rst$/,
		:zmat => /.*\.g?zmat$/,
		:sdf => /.*\.sdf$/,
		:mol2 => /.*\.mol2$/,
		:dftb_gen => /^in.gen$/
	}

	# File format modules
	require "classes/geometry_file_formats/geofile_pdb.rb"
	require "classes/geometry_file_formats/geofile_pdb.rb"
	require "classes/geometry_file_formats/geofile_xyz.rb"
	require "classes/geometry_file_formats/geofile_turbomole_coord.rb"
	require "classes/geometry_file_formats/geofile_amber_crd.rb"
	require "classes/geometry_file_formats/geofile_zmatrix.rb"
	require "classes/geometry_file_formats/geofile_sdf.rb"
	require "classes/geometry_file_formats/geofile_mol2.rb"
	require "classes/geometry_file_formats/geofile_dftb_gen"

	def GeoFile.read(geometry, filetype, arguments = {})
		#: Read file of specified filetype using the appropriate module
		Cuby::error("Unknown file format :#{filetype} in GeoFile.read") unless TYPE_MODULE_R.has_key?(filetype)
		typemodule = eval(TYPE_MODULE_R[filetype]) # Get module by evaluating its name
		typemodule.read(geometry, arguments)
	end

	def GeoFile.write(geometry, filetype = :xyz, arguments = {})
		#: Write file of specified filetype using the appropriate module
		Cuby::error("Unknown file format :#{filetype} in GeoFile.write") unless TYPE_MODULE_W.has_key?(filetype)
		typemodule = eval(TYPE_MODULE_W[filetype]) # Get module by evaluating its name
		typemodule.write(geometry, arguments)
	end

	def GeoFile.type_from_filename(filename)
		#: Guess file type from a file name (uses regular expressions defined in TYPE_FN_PATTERNS)
		#: Unknown file format returns nil.
		TYPE_FN_PATTERNS.each_pair{|key,pattern|
			return key if File.basename(filename).downcase =~ pattern
		}
		return nil
	end

	def GeoFile.type_from_string(string)
		#: Get file type (symbol) from its string representation (i.e. user input)
		TYPE_FN_PATTERNS.each_key{|key|
			return key if string.downcase == key.to_s
		}
		Cuby::error("Unknown geometry file type '#{string.downcase}'")
	end

	def  GeoFile.type_autodetect(filename)
		#: File type autodetection - reads file, detects format and returns pipe
		#: to the original contents
		rd, wr = IO.pipe
		f = File.open(filename,"r")
		filetype = nil
		counter = 0
		atomcount = 0
		xyzheader = false
		mol2 = false
		coord = false
		f.each {|line|
			atomcount += 1 if line =~ /^ATOM/
			xyzheader = true if counter == 0 && line =~ /^[0-9]+$/
			mol2 = true if line == "@<TRIPOS>MOLECULE"
			coord = true if counter == 0 && line =~ /^$coord/

			wr.print(line)
			counter += 1
		}
		if atomcount > 0
			filetype = :pdb
		elsif xyzheader
			filetype = :xyz
		elsif coord
			filetype = :tm_coord
		elsif mol2
			filetype = :mol2
		end

		wr.close
		f.close
		return [filetype, rd]
	end
end
