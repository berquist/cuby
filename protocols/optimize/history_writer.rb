
class HistoryWriter
	def initialize(geometry, history_fn, freq, selection = nil)
		# Save parameters
		@history_fn = history_fn
		@freq = freq
		@selection = selection
		if @selection
			@selection = geometry.atomlist_from_selection(@selection)
		end

		# No history written if freq == 0
		return if @freq == 0

		# Cycle counter
		@last_cycle_written = -1

		# Open the file
		@file = File.open(history_fn, "w+")
	end

	def write(cycle, geometry, values = {})
		return if @freq == 0
		# Write each @freq cycles
		if cycle % @freq == 0
			writestep(cycle, geometry, values)
			@last_cycle_written = cycle
		end
	end

	def write_final(cycle, geometry, values = {})
		return if @freq == 0
		# Write final geometry if it was not written yet
		unless @last_cycle_written == cycle
			writestep(cycle, geometry, values)
		end
	end

	def close
		@file.close if @file
	end

	#Private method that handles the actual writing
	def writestep(cycle, geometry, values = {})
		s = "Cycle: #{cycle}"
		if values[:energy]
			s += " E: #{values[:energy]} kcal/mol"
		end
		Cuby::log.puts_debug "Writing history to file #{@history_fn}"
		if @selection
			geometry.geometry_from_list(@selection).write_xyz(:file => @file, :second_line => s, :velocities => values[:velocities])
		else
			geometry.write_xyz(:file => @file, :second_line => s, :velocities => values[:velocities])
		end
		@file.flush
	end
	private :writestep
end

