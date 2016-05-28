require "classes/optimization/optimizer.rb"
require "classes/optimization/optimizer_common.rb"

class OptimizerTRIM< Optimizer
	include OptimizerCommon
	attr_reader	:hessian

	# Hessian estimate handling
	ALLOWED_H_ESTIMATE = [:number_inv, :vector, :matrix]
	
	def initialize(settings)
    		super(settings)

		@maxstep = @settings[:trust_radius]
		@maxstep_max = @settings[:trust_radius_max]
		@maxstep_min = @settings[:trust_radius_min]

		@failure_count = 0
	end

	def hessian_estimate_default
		# Default hessian initialization
		@hessian = Matrix.unit(@vector.size)
	end

	def initialize_hessian(argument)
		# Initialize hessian
		if argument.kind_of?(Numeric)
			# Size provided, use an universal diagonal guess
			@hessian = Matrix.unit(@vector.size) * argument
		elsif argument.kind_of?(Vector)
			# Vector provided, put it on diagonal
			@hessian = Matrix.diagonal(argument)
		elsif argument.kind_of?(Matrix)
			@hessian = argument
		else
			raise "Argument of initialize_hessian must be either a number, vector or a matrix"
		end
	end

	def newton_step_shift(g, h, shift)
		step = Vector.of_size(h.size)
		step.each_index{|i| step[i] = -1.0 * g[i] / (h[i] + shift)}
		return step
	end

	def step_vector_trim
		# Diagonalize the hessian
		evec, real, imag = @hessian.eigensystem

		# Diagonalized hessian
		dia_h = real.to_vector
		# Gradient in diagonal basis
		dia_g = evec.transpose * @gradient

		# Count eigenvalues
		thr = 1.0e-6
		c_zero = 0
		c_negative = 0
		nonzero_ev = []
		dia_h.each_index{|i|
			if dia_h[i] < -thr
				c_negative += 1
				nonzero_ev << dia_h[i]
			elsif dia_h[i] < thr
				c_zero += 1
				# Zero it
				dia_h[i] = 0.0
			else
				nonzero_ev << dia_h[i]
			end
		}

		# TS opimization options
		#!#

		lowest_ev = dia_h.min

		# Trust radius enforcement
		step = newton_step_shift(dia_g, dia_h, 0.0)
	
		if step.abs < @maxstep
			###puts_debug "Step size OK"
		else
			##puts_debug "Step size has to be scaled"
			# More-Sorensen
			shift = -lowest_ev + 0.0000001 #lowest_ev
			shift = 0.0 if shift < 0
			maxit = 100
			maxit.times{|i|
				raise "Trust radius solution failed, shift = #{shift}" if i == maxit - 1

				step = newton_step_shift(dia_g, dia_h, shift)
				stepsize = step.abs
				if (@maxstep - stepsize).abs < 0.00001
					#puts "Iterations to get shift: #{i}"
					break 
				end
				# shifted_hi = h - shift * I
				shifted_hi = Vector.of_size(dia_h.size)
				shifted_hi.each_index{|j| shifted_hi[j] = 1.0/(dia_h[j] + shift)}
				stepsize_deriv = -((shifted_hi.elementwise_multiply(step)) * step) / stepsize
				shift +=  stepsize / stepsize_deriv * (@maxstep - stepsize) / @maxstep
			}
		end

		step_in_orig = evec * step
		return step_in_orig
	end

	def step
		# Save current data
		@energy_p = @energy
		@grad_p = @gradient
		@hessian_p = @hessian

		# Call the calculation, printing and convergence check
		if @reuse_point
			@energy, @gradient = @reuse_point
			@reuse_point = nil
		else
			@calculation_queue_proc.call(0,@vector)
			@energy, @gradient = @calculation_run_proc.call(0)
			@stat_counters[:energy_gradient_calls] += 1
		end

		# Energy change
		@delta_e = @energy - @energy_p if @energy_p

		# Print step
		@print_step_proc.call

		# Check convergence
		if @cycle > 0
			return :converged if @check_convergence_proc.call(@delta_e) 
			@stat_counters[:steps_upward] += 1 if @delta_e > 0.0
		end

		# Trust radius update (from OptimizerCommon)
		step_status = tr_update
		return :failed if step_status == :failed

		# Hessian update (from OptimizerCommon)
		hessian_update_qn unless step_status == :refused

		# Propose step
		@step = fix_periodicity(step_vector_trim)

		# Predicted energy change
		@predicted_de = @gradient * @step + 0.5 * @hessian.inverse.vt_self_v(@step,@step)

		# Apply step
		@vector_p = @vector
		@vector = fix_periodicity(@vector + @step)

		# Success
		return :continue
	end
end
