require "classes/optimization/optimizer.rb"

#===============================================================================
# Nelder-Mead simplex optimization - gradient-free
#===============================================================================

class OptimizerNM < Optimizer
	INIT_DISPLACEMENT = 0.1

	class NMPoint
		attr_reader :coordinates
		attr_reader :energy

		def initialize(coordinates, energy)
			@coordinates = coordinates
			@energy = energy
		end
	end

	def initialize(settings)
		# Generic initialization
		super(settings)

		# Array of points to be stored
		@points = []

		# Parameters
		@alpha = 1.0
		@gamma = 2.0
		@rho = -0.5
		@sigma = 0.5
	end

	def hessian_estimate_default
	end
	
	def initialize_hessian(argument)
		# Do nothing...
	end

	def what_to_calculate
		# By default, optimizer needs energy and gradient.
		# Here, only energy is needed
		# ancestor optimizer.
		return [:energy]
	end
	
	def step
		if @cycle <= @vector.size 
			# In first N+1 points, build the initial simplex
			# Start with starting point + points displaced by 10%
			if @cycle == 0
				# First point in initial coordinates
				@starting_vector = @vector
			else
				displaced_c = @cycle - 1
				@vector = @starting_vector * 1
				@vector[displaced_c] = @vector[displaced_c] * (1.0 - INIT_DISPLACEMENT)
			end

			@calculation_queue_proc.call(0,@vector)
			@energy = @calculation_run_proc.call(0)[0]
			@points << NMPoint.new(@vector, @energy)

			@print_step_proc.call
			Cuby::log.puts_debug "NM optimizer: Filling points array"
			return :continue
		else
			# Order existing points so that first has the lowest energy
			@points.sort!{|a,b| a.energy <=> b.energy}
			best_e = @points[0].energy

			# Centroid of the simplex without the last point
			x0 = Vector.zero(@points[0].coordinates.size)
			(@points.size - 1).times.each{|i|
				x0 += @points[i].coordinates / (@points.size - 1)
			}

			# Index of last (n+1) point
			last_i = @vector.size

			# Calculate the reflected point
			xr = x0 + @alpha*(x0 - @points.last.coordinates)
			@calculation_queue_proc.call(0,xr)
			fxr = @calculation_run_proc.call(0)[0]
			Cuby::log.puts_debug "NM optimizer: Reflected point calculated"

			# Decision tree
			if @points[0].energy <= fxr && fxr < @points[last_i - 1].energy
				# reflected point is better than the second worst,
				# but not better than the best
				@points[last_i] = NMPoint.new(xr, fxr)
			elsif  fxr < @points[0].energy
				# The new point is the best so far
				# Calculate expanded point
				xe = x0 + @gamma * (x0 - @points.last.coordinates)
				@calculation_queue_proc.call(0,xe)
				fxe = @calculation_run_proc.call(0)[0]
				Cuby::log.puts_debug "NM optimizer: Expanded point calculated"
				if fxe < fxr
					# Expanded point is better than the reflected point
					@points[last_i] = NMPoint.new(xe, fxe)
				else
					@points[last_i] = NMPoint.new(xr, fxr)
				end
			else
				# Contraction
				Cuby::log.puts_debug "NM optimizer: Contracted point calculated"
				xc = x0 + @rho * (x0 - @points.last.coordinates)
				@calculation_queue_proc.call(0,xc)
				fxc = @calculation_run_proc.call(0)[0]
				if fxc < @points[last_i].energy
					@points[last_i] = NMPoint.new(xc, fxc)
				else
					# Reduction
					# Replace all the points but the best
					Cuby::log.puts_debug "NM optimizer: Reduction - recalculating most points"
					(1..last_i).each{|i|
						xi = @points[0].coordinates + @sigma  * (@points[i].coordinates - @points[0].coordinates)
						@calculation_queue_proc.call(0,xi)
						fxi = @calculation_run_proc.call(0)[0]
						@points[i] = NMPoint.new(xi, fxi)
					}
				end
			end

			# Sort the points again
			# Store the best of them as the global results
			@points.sort!{|a,b| a.energy <=> b.energy}
			@energy = @points[0].energy
			@vector = @points[0].coordinates
			@print_step_proc.call

			# Check convergence
			if @cycle > @vector.size
				@delta_e = @energy - best_e
				if @delta_e != 0
					Cuby::log.puts_debug "NM optimizer: DeltaE: #{@delta_e}"
					return :converged if @check_convergence_proc.call(@delta_e) 
				end
			end
			
			# Success
			return :continue
		end
	end
	
end
