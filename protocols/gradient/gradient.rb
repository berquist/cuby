
module ProtocolGradient
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :elementary_calculation
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Gradient calculation')
	end

	def prepare_protocol
		# Read geometry
		load_geometry
		# Create calculation
		@calculation = Calculation.new(@name, @settings, @geometry)
		# Prepare calculation
		what = [:energy, :gradient]
		what << :point_charges_gradient if @settings[:gradient_on_point_charges]
		@calculation.prepare(what)
	end

	def run_protocol(queue)
		# Submit calculation
		@calculation.send_to_queue(queue)
	end

	def print_protocol
		# Print results - energy
		results.print_energy
		Cuby::log.puts_v_only(:normal,"")
		results.print_energy_decomposition if @settings[:print].include?(:energy_decomposition)
		# Print results - gradient
		Cuby::log.puts_v_only(:brief, "Gradient:")
		Cuby::log.puts_v_only(:normal, "\nGradient (kcal/mol/A):")
		Cuby::log.puts(results.gradient.to_s)
		# Print results - gradient RMS
		Cuby::log.puts_v_only(:normal,"")
		Cuby::log.puts_v(:brief, "RMS grad: #{'%.6f' % results.gradient.rms} kcal/mol/A")
		if results.gradient_components.size > 0
			Cuby::log.puts
			Cuby::log.puts"Gradient decomposition:"
		end
		results.gradient_components.each_pair{|name, grad|
			Cuby::log.puts(name.to_s)
			Cuby::log.puts(grad.to_s)
		}

		# Gradient of atomic charges
		if results.atomic_charges_gradient
			Cuby::log.puts
			Cuby::log.puts "Gradient of atomic charges (elementary_charge/A):"
			Cuby::log.puts(results.atomic_charges_gradient.to_s)

		end

		# Gradient on point charges
		if results.point_charges_gradient
			Cuby::log.puts
			Cuby::log.puts "Gradient on point charges (kcal/mol/A):"
			Cuby::log.puts(results.point_charges_gradient.to_s)
		end

		# Optional output file - for gaussian as a driver
		if @settings.set?(:gaussian_external_output)
			require "classes/data_file_formats/gaussian_external_output.rb"
			GaussianExternalOutput.write(@settings[:gaussian_external_output], results, [:energy, :gradient]) 
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
			@results.gradient.project_out_trans_rot!(@geometry) if @settings[:gradient_projected]
		end
		return @results
	end
end
