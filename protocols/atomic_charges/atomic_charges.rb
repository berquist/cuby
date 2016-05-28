
module ProtocolAtomicCharges
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :elementary_calculation
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:fragment_N,          :optional, "Blocks defining fragments, N=(1..X)"]
	]
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Atomic charges calculation')
	end

	def prepare_protocol
		# Read geometry
		load_geometry
		# Create calculation
		@calculation = Calculation.new(@name, @settings, @geometry)
		# Prepare calculation
		@calculation.prepare([:energy, :atomic_charges])
	end

	def run_protocol(queue)
		# Submit calculation
		@calculation.send_to_queue(queue)
	end

	def print_protocol
		# Print results
		results.print_energy
		Cuby::log.puts_v_only(:normal,"")
		results.print_energy_decomposition if @settings[:print].include?(:energy_decomposition)
		results.print_atomic_charges(@geometry)
		results.print_dipole if @settings[:print].include?(:dipole)

		# Print charges on fragments, if fragments are present
		print_fragment_charges

		# Save the atomic charges
		if @settings.set?(:atomic_charges_write)
			results.atomic_charges.write_file(@settings[:atomic_charges_write])
		end
		# Save as point charges
		if @settings.set?(:atomic_charges_write_pch)
			results.atomic_charges.write_file_xyzc(@settings[:atomic_charges_write_pch], @geometry)
		end
	end

	def cleanup_protocol
		@calculation.cleanup
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			# Get results
			@results = @calculation.results
		end
		return @results
	end

	#=======================================================================
	# Private
	#=======================================================================
	def print_fragment_charges
		return unless @settings.has_block?(:fragment_1)
		Cuby::log.puts_v_only(:normal,"")
		Cuby::log.puts_v_only(:normal,"Charge on fragments:")
		i = 1
		while @settings.has_block?("fragment_#{i}".to_sym)
			sum = 0.0
			atomlist = @geometry.atomlist_from_selection(@settings["fragment_#{i}".to_sym, :selection])
			atomlist.each{|chi|
				sum += @results.atomic_charges[chi]
			}
			Cuby::log.puts_v_only(:normal, sprintf("%4d%10.4f", i, sum))
			i += 1
		end
	end
end
