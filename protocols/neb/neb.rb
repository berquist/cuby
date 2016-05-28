require "protocols/optimize/translation_rotation_removal.rb"
require "classes/geometry/geometry_series.rb"

module ProtocolNeb
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :driver
	PROTOCOL_DESCRIPTION = "Nudged elastic band method for minimum energy path optimization"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation, :required, "Setup for the calculation"]
	]
	#=======================================================================


	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Nudged Elastic Band')
	end

	def prepare_protocol
		# Load geometries - series in .xyz file needed for now
		unless FileTest.file?(File.expand_path(@settings[:geometry]))
			Cuby::error("For the NEB job, geometry keyword must reference to a .xyz file")
		end
		@geometries = GeometrySeries.new(@settings)
		if @geometries.size < 3
			Cuby::error("A series of at least three geometries is needed for the NEB job")
		end

		# Optimizer setup
		require "classes/optimization/optimizer_LBFGS.rb"
		@optimizer = OptimizerLBFGS.new(@settings)
		
		# Disable energy criterion in convergence check
		@settings[:opt_convlimit_e] = 0.0


		# Create calculations
		@calculations = []
		@geometries.each_index{|i|
			@calculations[i] = Calculation.new(@name + "_#{i}", @settings.block(:calculation), @geometries[i])
			@calculations[i].prepare([:energy, :gradient])
		}
	end

	def run_protocol(queue)
		# Translation/rotation remover
		transrot = TranslationRotationRemoval.new(@settings.block(:calculation))

		# Starting vector
		@optimizer.starting_vector = geometries_to_vector

		# Initialization of the hesian
		@optimizer.initialize_hessian(@settings[:opt_diagonal_h0]) if @optimizer.respond_to?(:initialize_hessian)

		# Calculation procedure: queue
		@optimizer.calculation_queue_proc = Proc.new{|calc_id, vector|
			puts vector.class
			$stdout.flush
			# Run the calculation (calc_id is ignored, assume serial optimizer only)
			update_geometries_from_vector(vector)
			# Queue calculations
			@calculations.each{|calc|
				calc.send_to_queue(queue)
			}
			# The actual vector on which the calculation ran must be returned
                        # Here it is identical to the input vector
			vector
		}

		# Calculation procedure: run
		@optimizer.calculation_run_proc = Proc.new{|calc_id|
			# Remove translation/rotation from gradient
			@calculations.each_index{|i|
				transrot.fix_gradient!(@geometries[i], @calculations[i].results.gradient)
			}

			# Composite energy - sum of energies
			#!# Does not include the elastic band
			energy = 0.0
			@calculations.each_index{|i|
				energy += @calculations[i].results.energy
			}

			# gradient = composite_gradient_vector
			gradient = elastic_band_gradient

			# Pass results to the optimizer
			[ energy, gradient ]
		}
		
		@optimizer.print_step_proc = Proc.new{
			if @settings[:optimize_print].include?(:steps)
				Cuby::log.puts "--------------------------------------------------------------------------------"
				Cuby::log.puts "Cycle: #{@optimizer.cycle}"
				Cuby::log.puts "Energies: (kcal/mol)"
				@calculations.each{|calc| Cuby::log.puts(sprintf("   %.4f", calc.results.energy))}

				# Write geometries
				write_path
			end
		}

		@optimizer.check_convergence_proc = Proc.new{|delta_e|
			convergence_check(delta_e, composite_gradient_vector)
		}


		# Run the optimizer
		@opt_status = @optimizer.run

	end

	def results
		# Return results of the last calculation
		return @calculations[0].results
	end

	def print_protocol
		case @opt_status
		when :converged
			puts "Converged in #{@optimizer.cycle + 1} steps"
		else
			puts "Optimizer finsihed with status: #{@opt_status}"
		end

		Cuby::log.puts ""
	end

	def cleanup_protocol
		@calculations.each{|calculation|
			calculation.cleanup
		}
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def convergence_check(delta_e, gradient)
		Cuby::log.puts "Convergence:" if @settings[:optimize_print].include?(:steps)

		# Convergence critera
		criteria = []
		criteria << {:name => "rmsGrad:", :value => gradient.rms, :limit => @settings[:opt_convlimit_rms_g] * @settings[:opt_quality]} if @settings[:opt_convlimit_rms_g] != 0.0
		criteria << {:name => "maxGrad:", :value => gradient.max_abs, :limit => @settings[:opt_convlimit_max_g] * @settings[:opt_quality]} if @settings[:opt_convlimit_max_g] != 0.0

		converged = true
		count = 0;
		criteria.each{|c|
			cond = c[:value].abs < c[:limit]
			if @settings[:optimize_print].include?(:steps)
				Cuby::log.puts sprintf("%s  %8.4f   limit: %8.4f  %s\n", c[:name], c[:value], c[:limit], cond ? "YES" : "NO")
			end
			converged = converged && cond
			count += 1 if cond
		}
		return converged
	end

	def write_path
		f = File.open(@settings[:neb_mep_file], "w+")
		@geometries.each_index{|i|
			second_line = "E: #{@calculations[i].results.energy} kcal/mol"
			@geometries[i].write_xyz(:file => f, :second_line => second_line)
		}
		f.close
	end
	
	def geometries_to_vector
		# Coordinates from all geometries merged into single long vector
		a = []
		@geometries.each{|geo|
			geo.to_vector.each{|x|
				a << x
			}
		}
		return Vector.from_array(a)
	end

	def update_geometries_from_vector(vector)
		# Updates all geometries from the merged vector
		vecindex = 0
		@geometries.each{|geo|
			geo.each{|atom|
				atom.x = vector[vecindex]
				atom.y = vector[vecindex + 1]
				atom.z = vector[vecindex + 2]
				vecindex += 3
			}
		}
		return nil
	end

	def composite_gradient_vector
		# Gradients from all calculations collected into one vector
		a = []
		@calculations.each{|calc|
			calc.results.gradient.to_vector.each{|x|
				a << x
			}
		}
		return Vector.from_array(a)
	end

	def elastic_band_gradient
		# The complete NEB gradient - gradient of images + the elastic band

		# Force constant
		fconst = @settings[:neb_fconst]

		# Get vectors of geometries and gradients for all images
		vecsize = @geometries[0].size * 3
		geo_vec = []
		grad_vec = []
		energies = []
		@geometries.each_index{|i|
			geo_vec[i] = @geometries[i].to_vector
			grad_vec[i] = @calculations[i].results.gradient.to_vector
			energies[i] = @calculations[i].results.energy
		}

		# Complete NEB gradient
		total_grad_vec = Vector.zero(@geometries.size * vecsize)

		# Copy unmodified gradient on end images to the complete gradient
		if @settings[:neb_optimize_endpoints]
			[0, @geometries.size - 1].each{|i|
				vecsize.times {|j|
					total_grad_vec[i*vecsize + j] = grad_vec[i][j]
				}
			}
		end

		# There are two versions of the NEB with different approach to correcting the kinks
		# Separate code according to @settings[:neb_version]

		# Iterate over all but end images		
		(1 .. (@geometries.size - 2)).each{|i|
			# Nudged elastic band: projecting out perpendicular band force and parallel image gradient

			# Tangent in point i
			if @settings[:neb_version] == :original
				# In this version, the tangent is estimated using central differences

				# Equation 1 in http://dx.doi.org/10.1063/1.1323224
				#tangent = geo_vec[i+1] - geo_vec[i-1]
				#tangent = tangent / tangent.abs

				# Better version of above - Equation 2 in http://dx.doi.org/10.1063/1.1323224
				tangent1= geo_vec[i] - geo_vec[i-1]
				tangent1 = tangent1 / tangent1.abs
				tangent2 = geo_vec[i+1] - geo_vec[i]
				tangent2 = tangent2 / tangent2.abs
				tangent = tangent1 + tangent2
				tangent = tangent / tangent.abs
			elsif @settings[:neb_version] ==  :improved
				# In this version, the tangent is constructed using the adjacent image
				# with higher energy
				# Equations 8 and 9 in http://dx.doi.org/10.1063/1.1323224
				t_plus = geo_vec[i+1] - geo_vec[i]
				t_minus = geo_vec[i] - geo_vec[i-1]
				if energies[i+1] > energies[i] && energies[i] > energies[i-1]
					tangent = t_plus
				elsif energies[i+1] < energies[i] && energies[i] < energies[i-1]
					tangent = t_minus
				else
					# In other cases, use both adjacent images, weighted
					# This ensures smooth switching between the above
					# Equations 10 and 11 in http://dx.doi.org/10.1063/1.1323224
					vmax = [(energies[i+1] - energies[i]).abs, (energies[i-1] - energies[i]).abs].max
					vmin = [(energies[i+1] - energies[i]).abs, (energies[i-1] - energies[i]).abs].min
					if energies[i+1] > energies[i-1]
						tangent = t_plus * vmax - t_minus * vmin
					else
						tangent = t_plus * vmin + t_minus * vmax
					end
				end
				# Normalize
				tangent = tangent / tangent.abs
			end

			# New gradient on i, perprendicular component removed
			# Equation 4 http://dx.doi.org/10.1063/1.1323224
			v = grad_vec[i] - grad_vec[i].dot(tangent) * tangent
			# Paste it into the total gradient
			vecsize.times {|j|
				total_grad_vec[i*vecsize + j] = v[j]
			}


			if @settings[:neb_version] == :original
				# Equation 5 in http://dx.doi.org/10.1063/1.1323224
				# Plain elestic band: force on this image
				fb = fconst * (geo_vec[i+1] - geo_vec[i]) - fconst * (geo_vec[i] - geo_vec[i-1])
				# Separte parallel band force
				fb_new = fb.dot(tangent) * tangent

				# Correction of kinks
				cosf = ((geo_vec[i+1] - geo_vec[i]).dot(geo_vec[i] - geo_vec[i-1])) /
					(geo_vec[i+1] - geo_vec[i]).abs /
					(geo_vec[i] - geo_vec[i-1]).abs
				swf = 0.5 * (1+Math.cos(Math::PI * cosf))
				fb_new += swf * (fb - fb_new)
			elsif @settings[:neb_version] == :improved
				# Equation 12 in http://dx.doi.org/10.1063/1.1323224
				fb = fconst * (geo_vec[i+1] - geo_vec[i]).abs - fconst * (geo_vec[i] - geo_vec[i-1]).abs
				fb_new = fb * tangent
			end

			# Add band force to total gradient
			vecsize.times {|j|
				total_grad_vec[j + i*vecsize] -= fb_new[j]
			}
		}
		return total_grad_vec
	end
end
