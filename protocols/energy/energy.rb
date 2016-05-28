
module ProtocolEnergy
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :elementary_calculation
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Energy calculation')
	end

	def prepare_protocol
		# Read geometry
		load_geometry
		# Create calculation
		@calculation = Calculation.new(@name, @settings, @geometry)
		# Prepare calculation
		what = [:energy]
		what << :mos if @settings[:print].include?(:molecular_orbitals)
		@calculation.prepare(what)
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
		results.print_dipole if @settings[:print].include?(:dipole)
		results.print_polarizability if @settings[:print].include?(:polarizability)

		# Orbitals
		if @settings[:print].include?(:molecular_orbitals)
			results.print_mos
		end

		# Optional output file - for gaussian as a driver
		if @settings.set?(:gaussian_external_output)
			require "classes/data_file_formats/gaussian_external_output.rb"
			GaussianExternalOutput.write(@settings[:gaussian_external_output], results, [:energy]) 
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
	# Optional protocol functionality
	#=======================================================================
	
	def single_line_result
		return sprintf("Energy: %20.6f kcal/mol", results.energy)
	end

	def single_number_result
		return results.energy
	end

	def update_geometry(geometry)
		# Update geometry, invalidate the results on previous one
		@geometry.copy_coordinates_from(geometry)
		@results = nil
	end

end
