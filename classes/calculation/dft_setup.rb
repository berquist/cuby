require "yaml"

class DftSetup

	# List of functionals known to cuby
	# lower case, without "-" characters
	FUNCTIONALS = {
		# Local
		"hfs"		=> :local,		# Hartree-Fock Slater
		"xalpha"	=> :local,		# The famous old Slater Xa theory
		"lsd"		=> :local,		# Local spin density (VWN-5A form)
		"vwn5"		=> :local,		# Local spin density (VWN-5)
		"vwn3"		=> :local,		# Local spin density (VWN-3)
		"pwlda"		=> :local,		# Local spin density (PW-LDA)
		# GGA
		"bnull"		=> :gga,		# Becke '88 exchange, no corr.
		"bvwn"		=> :gga,		# Becke '88 exchane, VWN-5 corr.
		"bp"		=> :gga,		# Becke X-Perdew 86 correlation
		"pw91"		=> :gga,		# Perdew-Wang GGA-II '91 func.
		"mpwpw"		=> :gga,		# Modified PW with PW correlation
		"mpwlyp"	=> :gga,		# same with LYP correlation
		"blyp"		=> :gga,		# Becke X with LYP correlation
		"gp"		=> :gga,		# Gill '96 X, Perdew '86 corr.
		"glyp"		=> :gga,		# Gill '96 X with LYP correlation
		"pbe"		=> :gga,		# Perdew-Burke-Ernzerhof
		"revpbe"	=> :gga,		# Revised PBE (exchange scaling)
		"rpbe"		=> :gga,		# Revised PBE (functional form of X)
		"pwp"		=> :gga,		# PW91 exchange + P86 correlation
		"olyp"		=> :gga,		# the optimized exchange and LYP
		"opbe"		=> :gga,		# the optimized exchange and PBE
		"xlyp"		=> :gga,		# the Xu/Goddard exchange and LYP
		"b97d"		=> :gga,
		# Meta-GGA
		"tpss"		=> :meta_gga,		# the TPPS functional
		# Hybrid
		"b1lyp"		=> :hybrid,		# One parameter Hybrid of BLYP
		"b3lyp"		=> :hybrid,		# Three parameter Hybrid of BLYP
		"b1p"		=> :hybrid,		# Analogous with Perdew exchange
		"b3p"		=> :hybrid,		# Analogous with Perdew exchange
		"g1lyp"		=> :hybrid,		# 1 par. analog with Gill 96 X
		"g3lyp"		=> :hybrid,		# 3 par. analog with Gill 96 X
		"g1p"		=> :hybrid,		# similar with P correlation
		"g3p"		=> :hybrid,		# similar with P correlation
		"pbe0"		=> :hybrid,		# 1 parameter version of PBE
		"pwp1"		=> :hybrid,		# 1 parameter version of PWP
		"mpw1pw"	=> :hybrid,		# 1 parameter version of mPWPW
		"mpw1lyp"	=> :hybrid,		# 2 parameter version of mPWLYP
		"pw91_1"	=> :hybrid,		# 1 parameter version of PW91
		"o3lyp"		=> :hybrid,		# 3 parameter version of OLYP
		"x3lyp"		=> :hybrid,		# 3 parameter version of XLYP
		"pw6b95"	=> :hybrid,		# Hybrid functional by Truhlar
		"bhlyp"		=> :hybrid,		# Becke Half and Half
		# Meta-hybrid
		"tpssh"		=> :meta_hybrid,	# hybrid version of TPSS with 10% HF exchange
		"tpss0"		=> :meta_hybrid,	# hybrid version of TPSS with 25% HF exchange
		"m06"		=> :meta_hybrid,
		"m062x"		=> :meta_hybrid,
		# Double-hybrid
		"b2plyp"	=> :double_hybrid,	# Grimme's 2006 double hybrid
		"mpw2plyp"	=> :double_hybrid,	# Schwabe/Grimme improved double hybrid
		"pwpb95"	=> :double_hybrid,	# New Grimme double hybrid
		# Range-separated
		"camb3lyp"	=> :range_separated,	# CAM-B3LYP
		"lcblyp"	=> :range_separated	# LC-BLYP
	}

	def initialize(functionals_file, grids_file = nil)
		#: Reads yaml files with definitions for each interface
		File.open(functionals_file) {|f| @functionals = YAML::load(f)}
		File.open(grids_file) {|f| @grids = YAML::load(f)} if grids_file
	end

	def functional(settings)
		#: Returns the functional name in format used by the interfaced program
		func = functional_internal(settings)
		return settings[:functional_custom] if func == 'custom'
		if settings[:functional].downcase == "custom"
			return settings[:functional_custom]
		else
			# Check whether it is available in the external program
			unless @functionals.has_key?(func)
				Cuby::error("Selected DFT functional not available in the interfaced program")
			end
			return @functionals[func]
		end
	end

	def grid(settings)
		#: Returns the functional name in format used by the interfaced program
		if settings[:dft_grid] == :custom
			return settings[:dft_grid_custom]
		else
			# Easy: the allowed values are checked by input parser
			return @grids[settings[:dft_grid]]
		end
	end

	def functional_type(settings)
		func = functional_internal(settings)
		if func == 'custom'
			return settings[:functional_custom_type]
		end
		return FUNCTIONALS[func]
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def functional_internal(settings)
		#: Returns the functional name in format used by cuby
		if settings[:functional].downcase == "custom"
			return 'custom'
		else
			# Remove extra characters, downcase
			func = settings[:functional].strip.downcase.gsub("-","")
			# Check if cuby knows it
			unless FUNCTIONALS.has_key?(func)
				Cuby::error("Selected DFT functional not recognized by cuby\nif you are sure you want to use it, use the 'custom' option")
			end
			return func
		end
	end
end
