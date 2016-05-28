class Regexp
	# Add safe match that can handle bad characters
	#!# This should be handled better!
	
	def match_safe(str)
		begin	
			result = match(str)
		rescue
			result = nil
		end

		return result
	end
end

class OutputParser

	class OutputParserError < RuntimeError
	end

	class ParserRecordCommon
		ALLOWED_GET	= [:line, :line_words, :matching, :submatch, :match_array, :number]
		ALLOWED_TYPE	= [:as_is, :float]

		def convert_type(result, type)
			case type
			when :float
				if result.kind_of?(Array)
					result.map!{|x| x.to_f}
				else
					result = result.to_f
				end
			when :as_is
				# Do nothing
			else
				raise "Wrong type to convert to"
			end
			return result
		end

		def line_to_result(line, get, matchdata)
			case get
			when :line
				result = line
			when :line_words
				result = line.strip.split
			when :matching
				result = matchdata[0]
			when :submatch
				result = matchdata[1]
			when :match_array
				result = []
				(matchdata.size - 1).times{|i| result[i] = matchdata[i+1]}
			when :number
				#!# Number matching should be improved to handle more possible notations
				if match2 = /[-]?[0-9]+\.[0-9]+/.match_safe(line)
					result = match2[0]
				else
					raise(OutputParserError, "Output parser expected a number in matched line, none found (#{@name})")
				end
			end
			return result
		end
	end

	class Line < ParserRecordCommon

		POSSIBLE_OPTIONS = [:get, :type, :multi, :required]
		ALLOWED_MULTI	= [:single, :first, :last]
		attr_reader :value

		def initialize(name, pattern, options = {}, parser_name = nil)
			# Type check
			raise("Name of searched item must be Symbol") unless name.kind_of?(Symbol)
			raise("Search pattern must be Regexp") unless pattern.kind_of?(Regexp)

			@name = name
			@pattern = pattern
			@parser_name = parser_name

			# Default options: get a number, return it as a float, raise error when multiple instances are found in the file
			@options = {
				:get => :number,
				:type => :float,
				:multi => :single,
				:required => true
			}
			@options.merge!(options)

			#!# Check allowed
			if (@options.keys - POSSIBLE_OPTIONS).size > 0
				raise "Unknown options present: #{(@options.keys - POSSIBLE_OPTIONS).join(', ')}"
			end

			@value = nil
			@count = 0
		end

		def apply(line)
			if matchdata = @pattern.match_safe(line)
				# Get the value
				result = line_to_result(line,  @options[:get], matchdata)

				# Convert type
				result = convert_type(result, @options[:type])

				# Check multiple occurence
				case @options[:multi]
				when :first
					@value = result if @count == 0
				when :last
					@value = result
				when :single
					if @count == 0
						@value = result 
					else
						raise(OutputParserError, "Multiple occurances of #{@name} found but only one should be present")
					end
				end
				@count += 1
				return true
			else
				return false
			end
		end

		def finalize
			if @options[:required]
				Cuby::error("Record #{@name} not found in the output#{@parser_name ? " (#{@parser_name})" : ""}") unless @value
			end
		end
	end

	class ErrorMsg
		def initialize(pattern, message, parser_name = nil)
			raise("Search pattern must be Regexp") unless pattern.kind_of?(Regexp)
			@pattern = pattern
			@message = message
			@parser_name = parser_name
			
		end

		def apply(line)
			if matchdata = @pattern.match_safe(line)
				message = @message
				message.gsub!("%1", matchdata[1]) if matchdata[1]
				Cuby::error(message + (@parser_name ? " (#{@parser_name})" : ""))
			end
		end
	end

	class Block < ParserRecordCommon
		POSSIBLE_OPTIONS = [:get, :type, :count]

		def initialize(name, pattern_start, pattern_valid, pattern_end, options = {}, parser_name = nil)
			# Type check
			raise("Name of searched item must be Symbol") unless name.kind_of?(Symbol)
			raise("Search pattern must be Regexp") unless pattern_start.kind_of?(Regexp)
			raise("Search pattern must be Regexp") unless pattern_valid.kind_of?(Regexp)
			raise("Search pattern must be Regexp") unless pattern_end.kind_of?(Regexp)

			@name = name
			@pattern_start = pattern_start
			@pattern_valid = pattern_valid
			@pattern_end = pattern_end
			@parser_name = parser_name

			# Default options: get a number, return it as a float, raise error when multiple instances are found in the file
			@options = {
				:get => :line,
				:type => :as_is,
				:count => false # (do not check item count)
			}
			@options.merge!(options)

			@count = 0
			@read = :wait
		end

		def apply(line)
			if @pattern_start.match_safe(line) && @read == :wait
				@read = :read 
			end
			@read = :finished if @pattern_end.match_safe(line) && @read == :read
			if @read == :read
				if matchdata = @pattern_valid.match_safe(line)
					# Get the value
					result = line_to_result(line,  @options[:get], matchdata)

					# Convert type
					result = convert_type(result, @options[:type])

					# Call code processing the record
					yield(@name, result, @count)
					@count += 1
				end
			end
		end

		def finalize
			if @options[:count]
				Cuby::error("Block of #{@options[:count]} records expected in the output, #{@count} found#{@parser_name ? " (#{@parser_name})" : ""}") unless @count == @options[:count]
			end
		end
	end

	def initialize(file, parser_name = nil)
		@file = file
		@errors = []
		@patterns = {}
		@blocks = {}
		@parser_name = parser_name
	end

	def add_pattern(name, pattern, options = {})
		@patterns[name] = Line.new(name, pattern, options, @parser_name)
	end

	def add_error(pattern, message)
		@errors << ErrorMsg.new(pattern, message, @parser_name)
	end

	def add_block(name, pattern_start, pattern_valid, pattern_end, options = {})
		@blocks[name] = Block.new(name, pattern_start, pattern_valid, pattern_end, options, @parser_name)
	end

	def execute
		# Open the file
		close = false
		if @file.class == String
			close = true
			@file = File.new(@file,"rb")
		end

		unless @file.kind_of?(IO)
			raise(TypeError,"File must be IO object (is #{@file.class})")
		end


		@file.each_line{|line|
			# Scan for errors
			@errors.each{|err| err.apply(line)}
			# Scan for values
			@patterns.each_pair{|name, pattern|
				pattern.apply(line)
			}
			# Scan for blocks
			@blocks.each_pair{|bname, block|
				block.apply(line){|name, data, count| yield(name, data, count)}
			}

		}

		# Close the file
		@file.close if close

		# Call final checks
		@patterns.each_pair{|name, pattern| pattern.finalize}
		@blocks.each_pair{|name, pattern| pattern.finalize}

		@executed = true
	end

	def [](name)
		raise "Parser must be executed prior to reading results" unless @executed
		return @patterns[name].value
	end
end
