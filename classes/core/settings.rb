################################################################################
#
# Class Settings
#
# Author: Jan Rezac
# Date created: 2008-11-03
# License: Cuby license
# Description: Settings is the object that handles cuby input
# Status: Works, Documented, may need some additions in working with subsections
#
################################################################################

#: This is the data structure that stores all the keywords specifying the job.
#: The indexes are Symbols, when the settings are read from hash or file, string
#: indexes are converted to symbols (and lowercase).

#: There is a list of possible keywords stored in YAML file with adittional info
#: on the keyword, such as type and default value.

#: When the Settings structure is build from a Hash (that can be read from YAML
#: file), every keyword is checked against the keyword list. Unknown keyword
#: raises an error, so does keyword that has apparently wrong type.

#: When a keyword is looked up using the [] method and the keyword is missing,
#: action specified in the keyword database is performed. The options are:
#* To die wit an error message (keyword is necessary)
#* To use default value and print a warning (where the default option might be questionable)
#* To use default value silently

#: When default value is used, it is written into the settings as well. It means
#: each warning is printed only once, and the complete settings can be printed
#: out at the end of the job.

#: The settings can be organized in blocks, and the blocks may be nested as
#: needed.

require "yaml"
require "classes/cuby.rb"

class SettingsKeyword
	#: One record used of the keywords database - there it is indexed by the keyword name.
	attr_accessor	:default_value
	attr_accessor	:when_not_present	# values: 
						# :die, :warning, :default, :other
	attr_accessor	:convert_to		# values: 
						# :integer, float, :string, :symbol, :boolean, :hash, :hash_symbol_symbol, :array, :symbol_list
						# :symbol_list - array of symbols
						# :hash_symbol_symbol - both keys and values are converted to lowcase symbols
	attr_accessor   :convert_options	# values:
						# :simple_array_of_strings - for convert_to=:array, allows simple input of list of strings
						# :keys_to_symbol - for convert_to=:hash, converts key to downcase symbols
	attr_accessor	:allowed_values		# array of allowed symbols
	attr_accessor	:range			# allowed range of numerical value
						#!# Not used yet
	attr_accessor   :commandline		# Short form of setting this keyword from commandline

	attr_accessor	:interface		# For interface-specific keywords: name of the interface
	attr_accessor	:protocol		# For protocol-specific keywords: name of the protocol
	attr_accessor   :devel			# Flag for development keyword not visible in the documentation
end

class Settings
	# class variables
	# @@keywords - stores the keyword database
	# @@blocks - stores the database of allowed block names - evaluated as regular expressions

	# internal variables
	attr_accessor :keywords
	attr_accessor :blocks

	attr_accessor :input_file # Name of the input file, if settings were created from a file

	#=======================================================================
	# Class initialization
	#=======================================================================
	
	def Settings.load_keyword_list(filename)
		#: Load keyword definition file (yaml file)
		Cuby::error "Keyword list file (#{filename}) not found" unless FileTest.exist?(filename)
		f = File.new(filename,"r")
		y = YAML::load(f)
		f.close
		@@keywords = y["keywords"]
		@@blocks = y["blocks"]
		@@hints = y["hints"]
		return true
	end

	def Settings.load_keyword_list_interfaces(path)
		#: Load interface-specific keywords

		# Initialize list of borrowed keywords
		borrowed_keywords = []
		
		# Process all interfaces
		Calculation.available_interfaces.each_pair{|name, dir|
			if FileTest.exist?(dir + "/keywords.yaml")
				f = File.new(dir + "/keywords.yaml")
				y = YAML::load(f)
				f.close
				# Tag interface keywords with interface name
				y["keywords"].each_pair{|kw_name, kw| 
					kw.interface = [] unless kw.interface
					kw.interface << name
				}
				# merge keywords
				@@keywords.merge!(y["keywords"]){|key, oldval, newval| Cuby::error "Interface '#{name}' attempts to owerwite definition of the keyword '#{key}'"}
				# merge blocks
				if y["blocks"]
					y["blocks"].each{|block|
						if @@blocks.include?(block)
							Cuby::error "Interface '#{name}' attempts to owerwite definition of a block '#{block.to_s}'"
						else
							@@blocks << block
						end
					}
				end
				# merge hints
				if y["hints"]
					y["hints"].each_pair{|hint, str|
						if @@hints.has_key?(hint)
							Cuby::error "Interface '#{name}' attempts to owerwite definition of a hint '#{hint.to_s}'"
						else
							@@hints[hint] = str
						end
					}
				end
				# borrowed keywords
				if y["borrowed_keywords"]
					y["borrowed_keywords"].each{|kw|
						kw = kw.to_sym
						borrowed_keywords << [kw, name]
					}
				end
			end
		}

		# Process borrowed keywords
		borrowed_keywords.each{|bk|
			kw, name = bk
			unless @@keywords[kw]
				raise "Problem with 'borrowed keyword' #{kw} in interface #{name} - keyword does not exist"
			end
			unless @@keywords[kw].interface
				raise "Problem with 'borrowed keyword' #{kw} in interface #{name} - it is not interface keyword"
			end
			@@keywords[kw].interface << name
		}
	end

	def Settings.load_keyword_list_protocols(path)
		#: Load protocols-specific keywords
		Job.available_protocols.each_pair{|name, dir|
			if FileTest.exist?(dir + "/keywords.yaml")
				f = File.new(dir + "/keywords.yaml")
				y = YAML::load(f)
				f.close
				# Tag interface keywords with interface name
				y["keywords"].each_pair{|kw_name, kw| kw.protocol = name}
				# merge keywords
				@@keywords.merge!(y["keywords"]){|key, oldval, newval| Cuby::error "Protocol '#{name}' attempts to owerwite definition of the keyword '#{key}'"}
				# merge blocks
				if y["blocks"]
					y["blocks"].each{|block|
						if @@blocks.include?(block)
							Cuby::error "Protocol '#{name}' attempts to owerwite definition of a block '#{block.to_s}'"
						else
							@@blocks << block
						end
					}
				end
				# merge hints
				if y["hints"]
					y["hints"].each_pair{|hint, str|
						if @@hints.has_key?(hint)
							Cuby::error "Protocol '#{name}' attempts to owerwite definition of a hint '#{hint.to_s}'"
						else
							@@hints[hint] = str
						end
					}
				end
			end
		}
	end

	def Settings.load_keyword_defaults(filename)
		#: Load yaml file containing modified keyword defaults. These are
		#: used for local configuration of cuby, the default defaults are
		#: listed in the keyword definition file.
		Cuby::error "Keyword defaults file (#{filename}) not found" unless FileTest.exist?(filename)
		defaults = Settings.from_file(filename)
		defaults.each_keyword{|kw|
			if @@keywords[kw].when_not_present == :default
				@@keywords[kw].default_value = defaults[kw]
			else
				Cuby::error("Configuration file (#{filename}) attempts to modify keyword that does not have default value")
			end
		}
	end
	
	#=======================================================================
	# Initialize
	#=======================================================================

	def initialize(parent = nil)
		#: Creates empty Settings object, optionaly with pointer to its
		#: parent in the tree of blocks of the input.
		@parent = parent
		@keywords = {}
		@blocks = {}
		@input_file = nil
	end

	def Settings.from_commandline(argv, require_input_file = true)
		options = {}
		files = []

		# Build a lookup for commandline switches
		short_to_keyword = {}
		@@keywords.each_pair{|name, keyword|
			if keyword.commandline
				short_to_keyword[keyword.commandline] = name
			end
		}

		# save argv, it is destroyed on processing
		argv = argv.dup

		while a = argv.shift
			# Is it option?
			if a =~ /^--.*/
				# Long format option
				option = a.downcase.sub(/^--/,'')
				# Change/remove charcters used by mistake
				option.gsub!("-","_")
				option.gsub!(/:$/,"")
				# Convert to symbol
				option = option.to_sym
				# Check if keyword exists
				unless @@keywords.has_key?(option)
					Cuby::error("Keyword '#{option}' set from commandline does not exist")
				end
			elsif a =~ /^-[a-z,A-Z,0-9]/
				# Short format option
				option = short_to_keyword[a]
				unless option
					Cuby::error("Unknown commandline option #{a}")
				end
			else
				option = nil
			end

			if option
				# Process option
				value = argv.shift
				Cuby::error("Value of commandline option #{a} missing") unless value
				value = Settings.fix_value(value, @@keywords[option], option)
				options[option] = value
			else
				# Otherwise, it is a file
				files << a
			end
		end

		# Process files
		if files.size == 0
			if require_input_file
				Cuby::error("One input file must be provided")
			else
				# Create empty settings
				settings = Settings.new
			end
		else
			if files.size > 1
				Cuby::error("Only one input file should be provided")
			end
			# Load settings
			settings = Settings.from_file(files[0])
			settings.input_file = files[0]
		end

		# Insert options into settings
		options.each_pair{|key, value|
			settings[key] = value
		}


		return settings
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def inspect
		#: Returns string with Hash representation of the settings built
		#: by -->to_hash method.
		return self.to_hash.to_yaml
	end

	def to_hash
		#: Converts the Settings to Hash representation in format
		#: recognized by -->Settings.from_hash
		hash = {}
		@keywords.each_pair{|name, keyword|
			hash[name.to_s] = keyword
		}
		@blocks.each_pair{|name, settings|
			hash[name.to_s] = settings.to_hash
		}
		return hash
	end

	def copy_from!(another_settings, overwrite, exceptions = [])
		overwrite_orig = overwrite
		case overwrite
		when :keep
			overwrite = false
		when :overwrite
			overwrite = true
		else
			raise "Invalid value of argument overwrite, use :overwrite or :keep"
		end

		exceptions_expanded = []
		exceptions.each{|exc|
			if exc.kind_of? Symbol
				# Add it directly
				exceptions_expanded << exc
			elsif exc.kind_of? String
				# Convert string to symbol
				exceptions_expanded << to_sym
			elsif exc.kind_of? Regexp
				# Find all matching keywords and blocks
				another_settings.keywords.each_key{|name|
					exceptions_expanded << name if exc.match(name.to_s)
				}
				another_settings.blocks.each_key{|name|
					exceptions_expanded << name if exc.match(name.to_s)
				}
			else
				raise "Keywords to be exceptrd should be specified by a Symbol, String or Regexp"
			end
		}

		another_settings.keywords.each_pair{|name, keyword|
			if (!set?(name) || overwrite) && !(exceptions_expanded.include?(name))
				if keyword.respond_to? :deep_copy
					@keywords[name] = keyword.deep_copy # Prevent deep linking, use a copy
				else
					@keywords[name] = keyword
				end
			end
		}

		another_settings.blocks.each_pair{|name, settings|
			if (!set?(name) || overwrite) && !(exceptions_expanded.include?(name))
				# Create block if needed
				new_block(name) unless @blocks[name]
				# merge block
				@blocks[name].copy_from!(settings, overwrite_orig, exceptions)
			end
		}
	end

	def child_block_as_copy(blockname, exceptions)
		#: Creates new Settings object as a child block of parent Settings and copies the
		#: parent contents
		self.new_block(blockname)
		child = self.block(blockname)
		exceptions2 = exceptions | [blockname] # Do not copy self
		child.copy_from!(self, :overwrite, exceptions2)
		return child
	end

	def new_block_as_copy_of(oldblock_name, newblock_name, exceptions = [])
		child = new_block(newblock_name)
		child.copy_from!(self.block(oldblock_name), :overwrite, exceptions)
		return child
	end

	#=======================================================================
	# Access to Keyword objects
	#=======================================================================
	
	def Settings.keyword(keyword_name)
		return @@keywords[keyword_name]
	end

	def Settings.all_keywords
		return @@keywords
	end

	#=======================================================================
	# Access
	#=======================================================================
	
	def [](*args)
		#: Read from the settings. Allows any depth of blocks, from root level:
		#| settings[:keyword]
		#: deeper into the block tree:
		#| settings[:block, :block, ..., :keyword]
		#: Behaviour when keyword is not found is described above.
		
		raise(ArgumentError,"Settings.[] needs at least one argument - the keyword.") if args.size < 1
		args.each{|key| raise(ArgumentError,"Settings.[] require Symbols as arguments") if key.class != Symbol }
		level = self
		string = ''
		(args.size-1).times {|section|
			level = level.blocks[args[section]]
			string += "#{args[section].to_s}:"
			Cuby::error("Block #{string} not found in the input.") if level == nil
		}
		r = level.keywords[args.last]
		if r == nil
			if @@keywords[args.last] == nil
				Cuby::error("Unknown keyword \"#{args.last.to_s}\" requested by cuby. This is weird...") 
			end

			case @@keywords[args.last].when_not_present
			when :die
				if level.keywords.size == 0
					empty = "\n(this input block is empty !!!)"
				else
					empty = ""
				end
				if level.parent
					level_name = level.parent.blocks.key(level)
					Cuby::error("Keyword \"#{args.collect{|i| i.to_s}.join(":")}\" not found in the input (in block #{level_name})."+empty)
				else
					Cuby::error("Keyword \"#{args.collect{|i| i.to_s}.join(":")}\" not found in the input (at the root level)."+empty)
				end
			when :warning
				Cuby::warning("Keyword \"#{args.collect{|i| i.to_s}.join(":")}\" not found in the input, using default value \"#{@@keywords[args.last].default_value}\".")
				r = @@keywords[args.last].default_value
				self[*args] = r
			when :default
				r = @@keywords[args.last].default_value
				self[*args] = r
			when :other
				raise "The code should handle keywords with when_not_present = :other (keyword #{args.collect{|i| i.to_s}.join(":")})"
			else
				raise("Wrong when_not_present value in keywords.yaml (#{args.last} has value #{@@keywords[args.last].when_not_present})")
			end
		end

		return r
	end

	def exists?(*args) # => boolean
		#: checks if requested keyword is present in settings
		raise(ArgumentError,"Settings.[] needs at least one argument - the keyword.") if args.size < 1
		args.each{|key| raise(ArgumentError,"Settings.[] require Symbols as arguments") if key.class != Symbol }
		level = self
		(args.size-1).times {|section|
			level = level.blocks[args[section]]
			return false if level == nil
		}
		return true if level.keywords[args.last] != nil
		return false
	end

	def set?(*args)
		return exists?(*args)
	end

	def []=(*args)
		#: Write to the settings - syntax similar to -->[] .
		#: Blocks are created when do not exist.
		#! settings[:block, :keyword] = value
		#
		raise(ArgumentError,"Settings.[]= needs at least two arguments - the keyword and value.") if args.size < 2
		(args.size-1).times {|i| raise(ArgumentError,"Settings.[]= require Symbols as arguments") if args[i].class != Symbol }
		level = self
		string = ''
		(args.size-2).times {|i|
			block = args[i]
			string += "#{block.to_s}:"
			unless level.blocks[block]
				raise "Attempt to create a block not defined in list of recognized blocks (#{string})" unless allowed_block?(block)
				level.blocks[block] = Settings.new(level)
			end
			level = level.blocks[block]
		}
		name = args[args.size - 2]
		raise("Attempt to write to a keyword not defined in keyword list'#{name.to_s}'") unless @@keywords.has_key?(name)
		value = args[args.size - 1]
		level.keywords[name] = value
	end

	def each_keyword(*args)
		#: Iterate over keywords in the given block (default: root)
		args.each{|key| raise(ArgumentError,"Settings.[] require Symbols as arguments") if key.class != Symbol }
		level = self
		(args.size).times {|section|
			level = level.blocks[args[section]]
			Cuby::error("Block #{args[section].to_s} not found in the input.") if level == nil
		}
		level.keywords.each_key{|key|
			# Kewords only
			yield key if @@keywords.has_key?(key)
		}
		return nil
	end

	def block(blockname)
		return @blocks[blockname]
	end

	def has_block?(blockname)
		return @blocks.has_key?(blockname)
	end

	def new_block(blockname)
		@blocks[blockname] = Settings.new(self)
		return @blocks[blockname]
	end

	def keyword_count
		return @keywords.size
	end

	def block_count
		return @blocks.size
	end

	def parent
		return @parent
	end

	def this_block_name
		#: Return name of the current block or :root if there is no
		#: parent object.
		if @parent
			return @parent.blocks.key(self)
		else
			return :root
		end
	end

	#=======================================================================
	# Save/load
	#=======================================================================

	def Settings.from_hash(hash, parent)
		settings = Settings.new(parent)

		hash.each_pair{|key, value|
			Cuby::error("Wrong keyword type in input (#{key.to_s})") unless key.class == String
			key = key.downcase.to_sym

			# Is it keyword?
			if kw = @@keywords[key]
				settings.keywords[key] = Settings.fix_value(value, kw, key.to_s)
				next
			end

			# Is it block?
			if settings.allowed_block?(key)
				settings.blocks[key] = from_hash(value, settings)
				next
			end

			if @@hints[key]
				Cuby::error("Problem found in the input:\n" + @@hints[key])
			end

			Cuby::error("Unknown entry '#{key.to_s}' found in the input")
		}
		return settings
	end

	def Settings.from_file(filename)
		#: Loads settings from YAML file (where it should be stored as a Hash)
		#: using all the fixes in -->from_hash

		@input_file = filename

		# check if file is ok
		if FileTest.exist?(filename)
			if FileTest.directory?(filename)
				Cuby::error("Input file \"#{filename}\" is not a file, but directory")
			end
		else
			Cuby::error("Input file \"#{filename}\" does not exist")
		end

		# open the file
		File.open(filename) {|f|

			begin
				h = YAML::load(f)
			rescue #Psych::SyntaxError
				#!# Psych is not used by versions < 1.9.1
				#!# Some better way to catch the exception should be implemented
				if $!.message =~ /line 0 column 0/
					Cuby::error("Wrong format of the input file '#{filename}',it seems it is not in YAML format")
				else
					m = $!.message.gsub(/.*line/,"line")
					Cuby::error("Syntax error in the input file '#{filename}' on #{m},\nplease make sure your input conforms to the YAML format")
				end
			end
			Cuby::error("Wrong format of input file, it should be YAML representation of hash.") unless h.class == Hash
			return Settings.from_hash(h, nil)
		}
	end

	def save_to_file(filename)
		#: Save yaml tree to file
		f = File.open(filename, "w+")
		f.puts self.to_hash.to_yaml
		f.close
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def allowed_block?(blockname)
		@@blocks.each{|regexp|
			return true if regexp.match(blockname.to_s)
		}
		return false
	end
	
	#=======================================================================
	# Class methods
	#=======================================================================
	
	def Settings.fix_value(value, keyword, key_name)
		#: Converts the value to the right type if needed or raises error.
		#: Key_name is used only in the error message.
		type = keyword.convert_to
		case type
		when :integer
			# Accepts integer, float with no decimal part and valid strings
			case value
			when Fixnum
				retval = value
			when Float
				if (value - value.floor) != 0.0
					Cuby::error("Illegal value of keyword \"#{key_name}\", expected integer.")
				else
					retval = value.to_i
				end
			when String
				if value =~ /^[+,-]{0,1}[0-9]+$/
					retval = value.to_i
				else
					Cuby::error("Illegal value of keyword \"#{key_name}\", expected integer.")
				end
			else
				Cuby::error("Illegal value of keyword \"#{key_name}\", expected integer.")
			end
			### Range check on retval
			return retval
		when :float
			# Accepts integer, float and valid strings
			case value
			when Fixnum
				retval = value.to_f
			when Float
				retval = value
			when String
				if value =~ /^[+,-]{0,1}[0-9]*\.[0-9]+[e]{0,1}[-]{0,1}[0-9]*$/
					retval = value.to_f
				else
					Cuby::error("Illegal value of keyword \"#{key_name}\", expected floating point number.")
				end
			else
				Cuby::error("Illegal value of keyword \"#{key_name}\", expected floating point number.")
			end
			# Range check on retval
			return retval
		when :string
			# Returns the string directly
			case value
			when String
				return value.chomp
			when Symbol
				return ":" + value.to_s
			else
				return value.to_s
			end
		when :symbol
			# Accepts string or symbol, both is converted to lowercase symbol
			case value
			when Symbol
				retval = value.to_s.downcase.to_sym
			when String
				retval =  value.downcase.to_sym
			else
				Cuby::error("Illegal value of keyword \"#{key_name}\", expected string or symbol.")
			end
			# Check for allowed values value
			if keyword.allowed_values
				unless keyword.allowed_values.include?(retval)
					Cuby::error("Keyword \"#{key_name}\" value \"#{retval}\" not recognized - see manual\nfor allowed options (e.g. using command \"cuby4 kw #{key_name}\")")
				end
			end
			return retval
		when :boolean
			# Anything that mean logical value can be used
			if value.class == TrueClass || value.class == FalseClass
				return value
			elsif value.class == String
				case value.downcase
				when "on"; return true
				when "true"; return true
				when "yes"; return true
				when "off"; return false
				when "false"; return false
				when "no"; return false
				else
					Cuby::error("Illegal value of keyword \"#{key_name}\", expected boolean value.")
				end
			else
				Cuby::error("Illegal value of keyword \"#{key_name}\", expected boolean value.")
			end
		when :hash
			if value.class == Hash
				if keyword.convert_options == :keys_to_symbol
					newhash = {}
					value.each_pair{|key, value| newhash[key.downcase.to_sym] = value}
					return newhash
				else
					return value
				end
			else
				Cuby::error("Illegal type of keyword \"#{key_name}\", expected an hash.")
			end
		when :hash_symbol_symbol
			if value.class == Hash
				converted = {}
				value.each_pair{|k,v|
					converted[k.to_s.downcase.to_sym] = v.to_s.downcase.to_sym
				}
				return converted
			else
				Cuby::error("Illegal type of keyword \"#{key_name}\", expected an hash.")
			end
		when :symbol_list
			case value
			when Symbol
				retval = [ value.to_s.downcase.to_sym ]
			when String
				retval = value.split(/\s*,\s*/).map{|x| x.downcase.to_sym}
			when Array
				retval = value.map{|x| x.to_s.downcase.to_sym}
			else
				Cuby::error("Illegal value of keyword \"#{key_name}\", expected list of symbols.")
			end
			# Check for allowed values value
			if keyword.allowed_values
				retval.each{|val|
					unless keyword.allowed_values.include?(val)
						Cuby::error("Keyword \"#{key_name}\" value \"#{val}\" not recognized - see manual for allowed options")
					end
				}
			end
			return retval
		when :array
			if value.class == Array
				case keyword.convert_options
				when :simple_array_of_strings
					# Array of strings: check if the items are strings
					value.each{|item|
						unless item.class == String
							Cuby::error("Illegal type of keyword \"#{key_name}\", expected an array of strings.")
						end
					}
				end

				return value
			else
				case keyword.convert_options
				when :simple_array_of_strings
					if value.class == String
						if value =~ /["']/
							Cuby::warning("Input parser attempts to interpret keyword \"#{key_name}\"
as a comma-separated list of strings, but the presence of quotation
marks indicate that the conversion might not be unambiguous.
Please use yaml array syntax for the list.")
						end
						return value.split(/\s*,\s*/)
					else
						Cuby::error("Illegal type of keyword \"#{key_name}\", expected an array of strings.")
					end
				else
					Cuby::error("Illegal type of keyword \"#{key_name}\", expected an array.")
				end
			end
		else
			raise("Wrong \"convert_to\" value for key \"#{key_name}\" in keywords.yaml")
		end
	end

	#=======================================================================
	# Helpers
	#=======================================================================
	
	def elements_hash(keyword)
		hash = {}
		unless keyword.class == Array
			keyword = [keyword]
		end
		self[*keyword].each_pair{|element, value|
			begin
				e = Atom.element_from_string(element)
			rescue
				Cuby::error "Keyword '#{keyword.map{|x|x.to_s}.join(":")}':\n#{$!.message}"
			end
			hash[e] = value
		}

		#!# Type checking removed for now

		return hash
	end

	#=======================================================================
	# Configuration keyword checks
	#=======================================================================

	def check_existence(keyword)
		if self[keyword] == "N/A"
			Cuby::error("Configuration keyword '#{keyword}' not set")
		end
		fn = File.expand_path(self[keyword])
		unless FileTest.exist?(fn)
			Cuby::error("Configuration keyword '#{keyword}' - invalid path:\n#{self[keyword]}")
		end
	end

	def check_file(keyword)
		check_existence(keyword)
		fn = File.expand_path(self[keyword])
		unless FileTest.file?(fn)
			Cuby::error("Configuration keyword '#{keyword}' should be a ifile:\n#{self[keyword]}")
		end
	end

	def check_dir(keyword)
		check_existence(keyword)
		fn = File.expand_path(self[keyword])
		unless FileTest.directory?(fn)
			Cuby::error("Configuration keyword '#{keyword}' should be a directory:\n#{self[keyword]}")
		end
	end

	def check_exe(keyword)
		check_file(keyword)
		fn = File.expand_path(self[keyword])
		unless FileTest.executable?(fn)
			Cuby::error("Configuration keyword '#{keyword}' should be executable:\n#{self[keyword]}")
		end
	end

	#=======================================================================
	# Helper classes - Seq - a sequence similar to the seq command
	#=======================================================================
	
	class Seq
		def initialize(start, step, stop)
			@start = start
			@step = step
			@stop = stop

			Cuby::error "seq(start, step, stop): start must be < stop" unless  start <= stop
			Cuby::error "seq(start, step, stop): step must be positive" unless  step > 0
		end

		def Seq.from_string(string)
			if matchdata = /^\s*seq\s*\(\s*([-]?[0-9]+\.?[0-9]*)\s*[,;]\s*([-]?[0-9]+\.?[0-9]*)\s*[,;]\s*([-]?[0-9]+\.?[0-9]*)\s*\)\s*/.match(string.downcase)
				return Seq.new(matchdata[1].to_f, matchdata[2].to_f, matchdata[3].to_f)
			else
				Cuby::error "Incorrect definition of a sequence, the syntax is 'seq(start, step, stop)'\n it is #{string}"
			end
		end

		def iterate
			x = @start
			while x <= @stop + 1.0e-8
				yield x
				x += @step
			end
		end

		def to_array
			a = []
			iterate{|x| a << x}
			return a
		end
	end

end
