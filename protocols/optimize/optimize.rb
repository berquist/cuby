require "protocols/optimize/translation_rotation_removal.rb"
require "protocols/optimize/history_writer.rb"
require "classes/misc/step_timer.rb"
require "classes/misc/dot_writer.rb"

# TBD:
# * Redundant coordinates - iterative decrease of internal step until transformation is successful and cartesian step is small

module ProtocolOptimize
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :driver
	PROTOCOL_DESCRIPTION = "Geometry optimization"
	#=======================================================================
	
	AUTO_RESTART_FILE = "auto_restart.xyz"

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Geometry optimization')
	end

	def prepare_protocol
		@timer = StepTimer.new(@settings)

		# Read geometry
		load_geometry

		if @settings[:opt_auto_restart] && FileTest.exist?(AUTO_RESTART_FILE)
			Cuby::log.puts("Automatic restart file found, updating coordinates")
			rst = Geometry.from_file(AUTO_RESTART_FILE)
			if rst.size == @geometry.size
				@geometry.copy_coordinates_from(rst)
				rst = nil
			else
				Cuby::warning("Number of atoms in restart file is different, it won't be used!")
			end
		end
		
		# Coordinate transformation setup
		@timer.step("transformation initialization") {   #!# seems to be extremely demanding !!!
			case @settings[:coordinates]
			when :cartesian
				if @settings[:optimize_region] != "%all()"
					# Can not be combined with :opt_freeze_cartesian

					require "protocols/optimize/coordinate_transformation_freezer.rb"
					@transformation = CoordinateTransformationFreezer.new(@geometry, @settings)
					@g2 = @geometry.geometry_from_selection(@settings[:optimize_region])

					# Switch off translation/rotation
					if @settings[:remove_translation]
						@settings[:remove_translation] = false
						Cuby::log.puts_debug("Freezer used, translation removal disabled")
					end
					if @settings[:remove_rotation]
						@settings[:remove_rotation] = false
						Cuby::log.puts_debug("Freezer used, rotation removal disabled")
					end
				elsif @settings.set?(:opt_freeze_cartesian)
					require "protocols/optimize/coordinate_transformation_cartesian_freezer.rb"
					@transformation = CoordinateTransformationCartesianFreezer.new(@geometry, @settings)

					# Switch off translation/rotation
					if @settings[:remove_translation]
						@settings[:remove_translation] = false
						Cuby::log.puts_debug("Freezer used, translation removal disabled")
					end
					if @settings[:remove_rotation]
						@settings[:remove_rotation] = false
						Cuby::log.puts_debug("Freezer used, rotation removal disabled")
					end
				else
					require "protocols/optimize/coordinate_transformation_none.rb"
					@transformation = CoordinateTransformationNone.new(@geometry, @settings)
				end
			when :"z-matrix"
				require "protocols/optimize/coordinate_transformation_zmatrix.rb"
				unless @geometry.info[:z_matrix]
					Cuby::warning("Automatically generated z-matrix will be used for optimization,\natom order will change.")
					zmat = ZMatrix.from_geometry(@geometry)
					@geometry = zmat.to_geometry
					@geometry.info[:z_matrix] = zmat
				end
				@transformation = CoordinateTransformationZMatrix.new(@geometry, @settings)
			when :redundant
				if @settings[:optimize_region] and @settings[:optimize_region] != "%all()"
					require "protocols/optimize/coordinate_transformation_redundant_freezer.rb"
					@transformation = CoordinateTransformationRedundantFreezer.new(@geometry, @settings)
					@g2 = @geometry.geometry_from_selection(@settings[:optimize_region])
				else
					require "protocols/optimize/coordinate_transformation_redundant.rb"
					@transformation = CoordinateTransformationRedundant.new(@geometry, @settings)
				end
			end
		}

		# Auto selection of optimizers
		if @settings[:optimizer] == :auto
			if @geometry.size <= 100
				@settings[:optimizer] = :quasi_newton
			else
				@settings[:optimizer] = :lbfgs
			end
		end

		# Prepare optimizer
		case @settings[:optimizer]
		when :trim
			require "classes/optimization/optimizer_TRIM.rb"
			@optimizer = OptimizerTRIM.new(@settings)
		when :rfo
			require "classes/optimization/optimizer_RFO.rb"
			@optimizer = OptimizerRFO.new(@settings)
		when :quasi_newton
			require "classes/optimization/optimizer_QN.rb"
			@optimizer = OptimizerQN.new(@settings)
		when :lbfgs
			require "classes/optimization/optimizer_LBFGS.rb"
			@optimizer = OptimizerLBFGS.new(@settings)
		when :plbfgs
			unless @settings[:coordinates] == :cartesian
				Cuby::error("PLBFGS optimizer can only be used with Cartesian coordinates!")
			end
			require "classes/optimization/optimizer_PLBFGS.rb"
			@optimizer = OptimizerPLBFGS.new(@settings)
		end

		# Geometry processing for geometry-specific optimizers
		if @optimizer.respond_to?(:get_geometry)
			@timer.step("geometry processing"){
				@optimizer.get_geometry(@transformation.geometry_to_geometry(@geometry))
			}
		end

		# Create calculation
		if @optimizer.parallel
			settings_0 = @settings.child_block_as_copy("calculation_para_#{0}".to_sym, [/^calculation_para/])
			@calculation = Calculation.new(@name, settings_0, @geometry)
		else
			# use settings
			@calculation = Calculation.new(@name, @settings, @geometry)
		end
		# Prepare calculation
		@calculation.prepare(@optimizer.what_to_calculate)


		# Arrays of calculations and their geometries
		@calculations = []
		@calc_geometries = []

		# The first item in the array is the main geometry and the associated calculation,
		# both accessible separately
		@calc_geometries << @geometry
		@calculations << @calculation

		# Build and prepare additional calculations for parallel optimizers
		if @optimizer.parallel
			Cuby::log.puts_debug "Parallel optimizer, preparing #{@optimizer.parallel} calculations"
			(@optimizer.parallel - 1).times{|i|
				@calc_geometries[i+1] = @geometry.deep_copy # Only complete copies are thread-safe
				settings = @settings.child_block_as_copy("calculation_para_#{i+1}".to_sym, [/^calculation_para/])
				@calculations[i+1] = Calculation.new(@name + "_#{i+1}", settings, @calc_geometries[i+1])
				@calculations[i+1].prepare(@optimizer.what_to_calculate)
			}
		end

		# Printing steps as dots
		if @settings[:optimize_print].include?(:steps_as_dots)
			@dot_writer = DotWriter.new
		end
	end

	def run_protocol(queue)
		# Initialize translation/rotation removal filter using values
		transrot = TranslationRotationRemoval.new(@settings)

		# Create file name for history
		if @settings.set?(:history_file)
			history_fn = @settings[:history_file]
		else
			if @settings.input_file
				history_fn = "history_" + File.basename(@settings.input_file).gsub(/\.[^\.]*$/,'') + ".xyz"
			else
				history_fn = "history.xyz"
			end
		end

		# History writers
		history = HistoryWriter.new(@geometry, history_fn, @settings[:history_freq])
		history_sel = HistoryWriter.new(@geometry, "sel_" + history_fn, @settings[:history_freq], @settings[:history_selection]) if @settings.set?(:history_selection)

		# Optimizer initialization

		# Starting vector
		@timer.step("initial coordinates transformation"){
			@optimizer.starting_vector = @transformation.geometry_to_vec(@geometry)
		}
		@optimizer.periodicity = @transformation.periodicity

		# Initialization of the hesian
		@optimizer.hessian_estimate_proc = Proc.new{
			if @optimizer.respond_to?(:initialize_hessian)
		 		@optimizer.initialize_hessian(@transformation.initial_hessian(@optimizer.stat_counters))
			end
		}

		# Calculation procedure: queue
		@optimizer.calculation_queue_proc = Proc.new{|calc_id, vector|
			# Run the calculation
			@timer.step("to cartesian"){
				@real_vector = @transformation.vec_to_geometry!(@calc_geometries[calc_id], vector)
			}
			@calculations[calc_id].send_to_queue(queue)
			@real_vector
		}

		# Calculation procedure: run
		@optimizer.calculation_run_proc = Proc.new{|calc_id|
			# Remove translation/rotation from gradient
			@timer.step("calculation"){
				transrot.fix_gradient!(@calc_geometries[calc_id], @calculations[calc_id].results.gradient)
			}

			# Pass results to the optimizer
			[
				@calculations[calc_id].results.energy,
				@timer.step("gradient transformation"){
					@transformation.gradient_to_vec(@calculations[calc_id].results.gradient)
				}
			]
		}
		
		# Checks step length in cartesians, returns step (if not in cartesian)
		if @settings[:coordinates] == :redundant
			@optimizer.check_cart_proc = Proc.new{|vector, step, ratio|   # ratio is another scale factor
				max_cart = 1.0
				old_cart_vec = @geometry.to_vector
				@transformation.vec_to_geometry!(@geometry, vector + step)
				cart_step = @geometry.to_vector - old_cart_vec
				@geometry.update_from_vector!(old_cart_vec)
				if cart_step.max_abs > max_cart
					Cuby::log.puts_debug "Stepsize downscaled by #{ratio * max_cart/cart_step.max_abs}"
					ratio * max_cart/cart_step.max_abs * step   # ratio -- so that it should not converge to nonzero vector
				else
					Cuby::log.puts_debug "Stepsize left as is."
					step
				end
			}
		end
		
		@optimizer.print_step_proc = Proc.new{
			if @settings[:optimize_print].include?(:steps)
				Cuby::log.puts "--------------------------------------------------------------------------------"
				Cuby::log.puts "Cycle: #{@optimizer.cycle}"
				Cuby::log.puts "Energy: #{@optimizer.energy}"
				# Print energy components
				@calculation.results.print_energy_decomposition if @settings[:print].include?(:energy_decomposition)
				@calculation.results.print_dipole if @settings[:print].include?(:dipole)
				@calculation.results.print_gradient if @settings[:print].include?(:gradient)
			end

			# Write history
			history.write(@optimizer.cycle, @geometry, :energy => @calculation.results.energy)
			history_sel.write(@optimizer.cycle, @geometry, :energy => @calculation.results.energy) if @settings.set?(:history_selection)

			if @settings[:opt_auto_restart]
				@geometry.write_xyz(:file => AUTO_RESTART_FILE)
			end
		}

		@optimizer.check_convergence_proc = Proc.new{|delta_e|
			@timer.step("gradient backtransformation"){
				cartgrad =  @transformation.transformed_gradient(@calculations[0].results.gradient.to_vector)
				if Cuby.log.logs[0].verbosity == :debug
					@g2 = @geometry.geometry_from_selection(@settings[:optimize_region]) unless @g2
					gradsizes = Vector.of_size(@g2.size)
					@g2.each_index{|i|
						gradsizes[i] = [cartgrad[i*3].abs, cartgrad[i*3+1].abs, cartgrad[i*3+2].abs].max
					}
					max = gradsizes.max

					@g2.each_index{|i|
						@g2[i].properties[:pdb_occupancy] = gradsizes[i]/max
					}
					@g2.write_pdb(:file => "grad.pdb", :extra_columns => true)
				end
				convergence_check(delta_e, cartgrad)
			}
		}


		# Run the optimizer
		@opt_status = @optimizer.run

		# Close the history
		history.write_final(@optimizer.cycle, @geometry, :energy => @calculation.results.energy)
		history.close
		if @settings.set?(:history_selection)
			history_sel.write_final(@optimizer.cycle, @geometry, :energy => @calculation.results.energy)
			history_sel.close
		end

		# Write the optimized/last geometry
		if @settings.set?(:restart_file)
			unless @settings[:restart_file] == ""
				filename = @settings[:restart_file]
				filetype = GeoFile.type_from_filename(filename)
				if filetype == :xyz
					s = "Energy #{@calculation.results.energy} kcal/mol"
					s += " Converged!" if @opt_status == :converged
					@geometry.write_file(filename, filetype, {:second_line => s, :append => @settings[:restart_file_append]})
				else
					@geometry.write_file(filename, filetype, {
						:append => @settings[:restart_file_append],
						:extra_columns => @settings[:pdb_extra_columns]
					})
				end
			end
		else
			s = "Energy #{@calculation.results.energy} kcal/mol"
			case @opt_status
			when :converged
				s += " Converged!" if @opt_status == :converged
				@geometry.write_xyz(:file => "optimized.xyz", :second_line => s, :append => @settings[:restart_file_append])
			else
				@geometry.write_xyz(:file => "last.xyz", :second_line => s, :append => @settings[:restart_file_append])
			end
		end
	end

	def results
		# Return results of the last calculation
		return @calculation.results
	end

	def print_protocol
		# Finish dot writer
		if @settings[:optimize_print].include?(:steps_as_dots)
			@dot_writer.finish
		end

		# Print final status
		if @opt_status == :converged
			steps = " in #{@optimizer.cycle + 1} steps"
		else
			steps = ""
		end
		Cuby::log.puts "Status: #{@opt_status}" + steps
		if @opt_status == :failed
			Cuby::log.puts ""
			Cuby::log.puts "!" * 80
			Cuby::log.puts "! OPTIMIZATION FAILED"
			Cuby::log.puts "!" * 80
			Cuby::log.puts "!"
			Cuby::log.puts "! Optimizer can not make a reasonable step although the step size limit"
			Cuby::log.puts "! had been decreased to the minimum. Try restarting the optimization"
			Cuby::log.puts "! from the last geometry."
			Cuby::log.puts "!"
			Cuby::log.puts "! If the restart fails again, it usually indicates that the gradient"
			Cuby::log.puts "! is not a correct/accurate enough derivative of the energy."
			Cuby::log.puts "!" * 80
			Cuby::log.puts ""
		end
		if @settings[:optimize_print].include?(:final_energy)
			@calculation.results.print_energy
		end
		@optimizer.print_stats if @settings[:optimize_print].include?(:statistics)
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
		if @settings[:opt_convlimit_e] != 0.0
			criteria << {:name => "deltaE: ", :value => delta_e, :limit => @settings[:opt_convlimit_e] * @settings[:opt_quality]} 
		else
			# print deltaE even when it is not evaluated
			puts sprintf("%s  %8.4f", "deltaE: ", delta_e)
		end
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
		if @settings[:optimize_print].include?(:steps_as_dots)
			@dot_writer.print_dot(['.','o','O','*'][count])
		end
		return converged
	end
end
