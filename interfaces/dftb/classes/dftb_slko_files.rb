
# This class handles the selection of the parameter set and determines
# the naming convention used for the parameter files.

class DftbSlkoFiles
	# Path to the Slater-Koster files directory
	attr_reader :path

	# Elements available in the set
	attr_reader :elements

	# Attributes needed to construct the filename from elements
	attr_reader :suffix
	attr_reader :separator
	attr_reader :downcase

	def initialize(settings)
		# Determine the path and file naming convention

		# Check the path
		unless FileTest.directory?(settings[:dftb_slko_basepath])
			Cuby::error("Path to the SLKO files (keyword dftb_slko_basepath) is not set to an existing\ndirectory", self)
		end
		@path = settings[:dftb_slko_basepath] + "/" + settings[:dftb_slko_set]
		Cuby::log.puts_debug "SLKO directory: #{@path}"
		unless FileTest.directory?(@path)
			Cuby::error("Selected set of SLKO files (#{settings[:dftb_slko_set]}) not found in\n#{settings[:dftb_slko_basepath]}", self)
		end

		# Get list of files
		files = Dir.entries(@path).sort - [".", ".."]
		files.delete_if{|fn| fn =~ /readme/i }

		# Determine file naming convention
		if settings.set?(:dftb_slko_format)
			determine_slko_naming(settings[:dftb_slko_format])
		else
			Cuby::log.puts_debug "Attempting to determine file naming of SLKO files automagically, dftb_slko_format not set"
			Cuby::log.puts_debug "First file used as a template: #{files[0]}"
			determine_slko_naming(files[0])
		end

		# Filter only slko files
		files.delete_if{|fn| fn !~ /#{@suffix}/}

		# List of elements for which parameters exist
		@elements = []
		if @separator != ""
			# If there is a separator, it is easy
			files.each{|fn|
				@elements |= [Atom.element_from_string(fn.gsub(/#{@separator}.*/,""))]
			}
		else
			# Without a separator, get elements from homodimers
			names = files.map{|s| s.gsub(@suffix,'').downcase}
			# Find homodimers - 2 or 4 characters
			names.each{|name|
				case name.size
				when 2
					@elements << Atom.element_from_string(name[0]) if name[0] == name[1]
				when 4
					@elements << Atom.element_from_string(name[0..1]) if name[0] == name[2] && name[1] == name[3]
				end
			}

		end

		return nil
	end

	def check_geometry(geometry)
		# Look for missing parameters
		unknown = geometry.elements_in_system - @elements
		if unknown.size > 0
			Cuby::error("Missing DFTB parameters (Slater-Koster files) for elements #{unknown.join(", ")}", self)
		end
		return nil
	end

	def determine_slko_naming(template)
		# Determine the SLKO naming convention from a template (actual file name)
		@suffix = /\.[^.]+$/.match(template)[0]
		if matchdata = /^[a-z,A-Z]+([^a-z,A-Z])*[a-z,A-Z]+\..*/.match(template)
			@separator = matchdata[1]
			@separator = "" if @separator == nil
		else
			@separator == ""
		end
		@downcase = template =~ /^[a-z]/ ? true : false
		return nil
	end

	def filename(element1, element2)
		# get full filename (including the path) from elements
		element1 = element1.to_s
		element2 = element2.to_s
		if @downcase
			element1.downcase!
			element2.downcase!
		else
			element1.upcase!
			element2.upcase!
		end
		return @path + "/" + element1 + @separator + element2 + @suffix
	end

end
