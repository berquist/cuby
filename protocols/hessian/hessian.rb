
module ProtocolHessian
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :elementary_calculation
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Hessian calculation')
	end

	def prepare_protocol
		# Read geometry
		load_geometry
		# Create calculation
		@calculation = Calculation.new(@name, @settings, @geometry)
		# Prepare calculation
		@calculation.prepare([:energy, :gradient, :hessian])
	end

	def run_protocol(queue)
		# Submit calculation
		@calculation.send_to_queue(queue)
	end

	def print_protocol
		# Print results - energy
		results.print_energy
		results.print_energy_decomposition if @settings[:print].include?(:energy_decomposition)
		# Gradient (optional)
		if @settings[:print].include?(:gradient)
			Cuby::log.puts_v_only(:normal,"Gradient (kcal/mol/A):")
			Cuby::log.puts_v_only(:normal, results.gradient.to_s)
			Cuby::log.puts_v_only(:normal,"")
		end
		# Print results - gradient RMS
		Cuby::log.puts_v_only(:normal,"")
		Cuby::log.puts_v(:brief, "RMS grad: #{'%.6f' % results.gradient.rms} kcal/mol/A")
		# Print results - hessian
		Cuby::log.puts_v_only(:normal,"")
		Cuby::log.puts_v_only(:normal,"Hessian (kcal/mol/A^2):")
		Cuby::log.puts(results.hessian.to_s("%10.4f"))

		# Save the hessian to a YAML file
		unless @settings[:hessian_write] == ""
			File.open(@settings[:hessian_write],"w+"){|f|
				f.puts (results.hessian).to_yaml
			}
		end

		# Optional output file - for gaussian as a driver
		if @settings.set?(:gaussian_external_output)
			require "classes/data_file_formats/gaussian_external_output.rb"
			GaussianExternalOutput.write(@settings[:gaussian_external_output], results, [:energy, :gradient, :hessian]) 
		end
		
		if @settings[:coordinates] == :redundant
			require "classes/objects/hessian_internal_coordinates.rb"
			require "protocols/optimize/coordinate_transformation_redundant.rb"
			
			transformation = CoordinateTransformationRedundant.new(@geometry, @settings)
			int_hessian = transformation.hessian_to_matrix(results.hessian, results.gradient)
			
			hess_int_coords = HessianInternalCoordinates.from_molecule(@geometry, transformation, int_hessian)
			hess_int_coords.write_yaml("hess_diag_internal.yaml")
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

			# Project hessian - optional
			if @settings[:hessian_projected]
				projector = Hessian.transrot_removal_projector(@geometry)
				@results.hessian = projector.transpose * @results.hessian * projector

				# Project gradient as well
				if @results.gradient
					@results.gradient.project_out_trans_rot!(@geometry)
				end
			end
		end

		return @results
	end
	
end
