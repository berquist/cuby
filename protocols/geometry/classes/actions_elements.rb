module ActionsElements
	def element_list
		# No geometry is written
		@write_geo = false

		# Get list
		list = @geometry.elements_in_system

		# Print it according to verbosity level
		Cuby::log.puts_v(:brief, object = 'Elements in the system:')
		Cuby::log.puts(list.join(", "))

	end

	def element_match
		# No geometry is written
		@write_geo = false

		# Get list
		in_geometry = @geometry.elements_in_system

		# Input list
		list = []
		@settings[:element_match_list].each{|e|
			e = e.to_s.downcase.capitalize.to_sym
			if PeriodicTable::ELEMENTS.include?(e)
				list << e
			else
				Cuby::error "Unknown element found on the list: #{e}"
			end
		}
		list.uniq!

		case @settings[:element_match_condition]
		when :at_least
			# matches molecules containing all listed elements and possibly others
			condition_matches = (in_geometry & list).size == list.size
		when :exactly
			#  matches molecules consisting of all listed elements but not more
			condition_matches = in_geometry.sort == list.sort
		when :no_more_than
			# matches molecules that contain some of the listed elements but not these outside the list
			condition_matches = ((in_geometry | list).size == list.size) && (in_geometry & list).size > 0
		when :some_from
			# matches molecules containg at lest one element from the list
			condition_matches = (in_geometry & list).size > 0
		end

		Cuby::log.puts_v(:brief, object = 'Elements in the system match the list:')
		Cuby::log.puts condition_matches ? "YES" : "NO"
	end
end
