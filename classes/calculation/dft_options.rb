class DftFunctional
	# Container for description of a DFT functional

	attr_reader :name # The name of the functional as recognized by the program
	attr_reader :type # Type - see FUNCTIONAL_TYPES
	attr_reader :description

	FUNCTIONAL_TYPES = [
		:lda,			# Local density approximation
		:gga,			# Generalized gradient approximation
		:meta_gga,		# Meta-GGA
		:hybrid,		# Hybrid GGA
		:meta_hybrid,		# Hybrid meta-GGA
		:double_hybrid, 	# Double-hybrid (PT2 correlation)
		:range_separated	# Long-range corrected
	]

	def initialize(name, type, description = "")
		# Value check
		raise("Unknown type of functional") unless FUNCTIONAL_TYPES.include?(type)
		# Initialize
		@name = name
		@type = type
		@description = description
	end
end

class DftOptions
	attr_reader :functionals
	# List of functionals. Index is name of the fucntional, 
	# in lower case, alphanuperic characters only

	attr_reader :grids

	def initialize
		@functionals = {}
		@grids = nil
	end

	#=======================================================================
	# Methods for accessing the options
	#=======================================================================
	
	def functional(settings)
		#: Returns the functional name in format used by the interfaced program
		return functional_data(settings).name

	end

	def functional_data(settings)
		func = settings[:functional].downcase.gsub(/[\+\-\(\)]/,'')

		# If it is a custom functional, add it to the list 
		# (so that other properties can be retrieved later)
		if func == "custom"
			Cuby::log.puts_debug "Adding custom DFT functional"
			@functionals['custom'] = DftFunctional.new(
				settings[:functional_custom],
				settings[:functional_custom_type], 
				"Custom functional")
		end

		# Return functional name
		if @functionals[func]
			return @functionals[func]
		else
			Cuby::error("DFT functional '#{settings[:functional]}' not available in interface #{settings[:interface]}")
		end
	end
end





