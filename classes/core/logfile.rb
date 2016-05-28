

class Logfile
	VERBOSITY_LEVELS = [
		:silent,
		:minimal,
		:brief,
		:normal,
		:debug
	]

	DEBUG_LINE_START = "#> "

	attr_reader :verbosity

	def initialize(verbosity_level = :normal, file = :stdout, errorfile = :stderr)
		@verbosity = verbosity_level
		@file_close = false
		to_file(file, errorfile)
	end

	def verbosity=(verbosity_level)
		if VERBOSITY_LEVELS.include?(verbosity_level)
			@verbosity = verbosity_level
		else
			raise("Unknown verbosity level set")
		end
	end

	def to_file(file, errorfile = :same)
		# Close old file when needed
		if @file_close
			@file_close = false
			@file.close
		end

		@file = nil

		# standard output
		if file == :stdout
			@file = $stdout
		end

		# Output to file specified by name
		if file.class == String
			@file = File.new(file,"w+")
			@file_close = true
		end

		# Output to IO object
		if file.kind_of?(IO)
			@file = file
		end

		# Check whether some output has been selected
		raise(TypeError, "Logfile not specified correctly (is #{file.class})") unless @file

		# Error file
		case errorfile
		when :stderr
			@file_err = $stderr
		when :same
			@file_err = file
		else
			raise(ArgumentError, "Logfile error output specification should be either :stderr or :same")
		end

		return true
	end

	def puts(object = '')
		@file.puts(object) unless @verbosity == :silent
	end	

	def pute(object = '')
		@file_err.puts(object)
	end	

	def flush
		@file.flush
		@file_err.flush
	end

	def puts_v(verbosity, object = '')
		@file.puts(object) if VERBOSITY_LEVELS.index(verbosity) <= VERBOSITY_LEVELS.index(@verbosity)
	end

	def puts_v_only(verbosity, object = '')
		if verbosity.kind_of?(Symbol)
			@file.puts(object) if verbosity == @verbosity || (verbosity == :normal && @verbosity == :debug)
		elsif verbosity.kind_of?(Array)
			@file.puts(object) if verbosity.include?(@verbosity) || (verbosity.include?(:normal) && @verbosity == :debug)
		else
			raise(TypeError,  "verbosity must be a symbol or an array")
		end
	end

	def puts_debug(object = '')
		@file.puts DEBUG_LINE_START + object.to_s if @verbosity == :debug
	end

	#=======================================================================
	# ASCII art
	#=======================================================================
	
	def print_cuby_logo(options = {})
		#: Prints Cuby logo in ascii art (when verbosity > verbosity_limit)
		#| print_cuby_logo('text')
		#|       _______
		#|      /\______\
		#|     / /      /
		#|    / / Cuby /   text
		#|    \/______/
		#|
		default_options = {
			:text => '',
			:verbosity_limit => :normal,
			:comment => false
		}
		options = default_options.merge(options) # Overwrite with defaults if no options present

		cs = ''
		cs = '#' if options[:comment]

		if VERBOSITY_LEVELS.index(options[:verbosity_limit]) <= VERBOSITY_LEVELS.index(@verbosity)
			@file.puts cs + '        _______  '
			@file.puts cs + '       /\______\ '
			@file.puts cs + '      / /      / '
			@file.puts cs + '     / / Cuby /  ' + ' ' + options[:text]
			@file.puts cs + '     \/______/   '
			@file.puts cs + '                 '
			return true
		end
		return false
	end

end
