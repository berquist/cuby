
class DefaultMethodSettings

	def initialize(lookup_tree)
		@lookup_tree = lookup_tree
	end

	def find(settings)
		# Finds and returns a Settings object containing defaults loaded
		# for a method specified in the Settings passed as an argument.
		# If none are found, returns false
	
		level = @lookup_tree
		kwlist = []
		Cuby::log.puts_debug "Searching for method-dependent default parameters:"

		while level["keyword"] || level["keyword_parent"]

			if level["keyword"]
				Cuby::log.puts_debug "   keyword #{level["keyword"]}"
				keyword = level["keyword"].to_s.downcase.to_sym
				if settings.set?(keyword)
					input_value = settings[keyword]
					Cuby::log.puts_debug "      found in input: #{input_value}"
				elsif [:default, :warning].include?(Settings.keyword(keyword).when_not_present)
					input_value = settings[keyword]
					Cuby::log.puts_debug "      default used: #{input_value}"
				else
					Cuby::log.puts_debug "      not set, aborting search"
					return false 
				end
				kwlist << [keyword.to_s, input_value]
			else # Check if there is a parent settings
				Cuby::log.puts_debug "   parent keyword #{level["keyword_parent"]}"
				return false unless settings.parent
				keyword = level["keyword_parent"].to_s.downcase.to_sym
				if settings.parent.set?(keyword)
					input_value = settings.parent[keyword]
					Cuby::log.puts_debug "      found in input: #{input_value}"
				elsif [:default, :warning].include?(Settings.keyword(keyword).when_not_present)
					input_value = settings.parent[keyword]
					Cuby::log.puts_debug "      default used: #{input_value}"
				else
					Cuby::log.puts_debug "      not set, aborting search"
					return false 
				end
				kwlist << ["parent." + keyword.to_s, input_value]
			end

			# Modify the value if requested
			if level["modify_string"]
				case level["modify_string"].downcase.to_sym
				when :downcase
					input_value.downcase!
				when :upcase
					input_value.upcase!
				end
			end
			if level["remove_substring"]
				input_value.gsub!(level["remove_substring"],"")
			end

			if level["values"].include?(input_value)
				level = level["values"][input_value]
			else
				# Value of the keyword not present in the parameter tree
				return false
			end
		end
		# level now points to the list of parameters

		# Create settings object from the hash
		return Settings.from_hash(level, nil)
	end

	def find_and_merge(settings)
		method_settings = find(settings)
		Cuby::log.puts_debug "Loading method-dependent default parameters:"
		if method_settings
		# Print debugging data
			# Copy defaults into the input when the keyword is not set
			method_settings.each_keyword{|kw|
				if !settings.set?(kw) && method_settings.set?(kw)
					settings[kw] = method_settings[kw]
					Cuby::log.puts_debug "   #{kw}: #{method_settings[kw]}"
				end
			}
			return true
		else
			Cuby::log.puts_debug "   none found"
			return false
		end
	end

end
