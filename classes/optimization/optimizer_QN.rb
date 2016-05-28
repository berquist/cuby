require "classes/optimization/optimizer.rb"
require "classes/optimization/optimizer_common.rb"

class OptimizerQN < Optimizer
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

		@step_scaling = @settings[:opt_tr_mode]
		@step_scaling = :components if @step_scaling == :default
	end
	
	def hessian_estimate_default
		# Default hessian initialization
		@hessian = Matrix.unit(@vector.size)
	end

	def initialize_hessian(argument)
		# Initialize hessian
		# in QN optimizer, inverse of the hessian is used
		if argument.kind_of?(Numeric)
			# Size provided, use an universal diagonal guess
			@hessian = Matrix.unit(@vector.size) * (1.0 / argument)
		elsif argument.kind_of?(Vector)
			# Vector provided, put it on diagonal
			@hessian = Matrix.diagonal(argument)
			# Invert the matrix
			@hessian.m.times{|i|
				@hessian[i,i] = 1.0 / @hessian[i,i]
			}
		elsif argument.kind_of?(Matrix)
			@hessian = argument.inverse
		else
			raise "Argument of initialize_hessian must be either a number, vector or a matrix"
		end
	end
	
	def step_vector_inverse_qn
		# Quasi-newton step using inverse hessian
		s = @hessian * -@gradient
		
		return s
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
		hessian_update_inverse_qn unless step_status == :refused

		# Propose step
		@step = fix_periodicity(step_vector_inverse_qn)

		# check step in cartesians
		if @check_cart_proc
			i = 0
			ratio = 0.9
			while @step != (s = fix_periodicity(@check_cart_proc.call(@vector, scale_step(@step, @maxstep), ratio)))
				@step = s
				i += 1
				ratio *= 0.7
				Cuby::log.puts_debug "#{i}. step of checking cartesian stepsize"
				break if i > 5
			end
		end

		@step = fix_periodicity(linesearch(@step))

		# Predicted energy change
		#!# The inverse should be replaced by another formula if possible!
		@predicted_de = @gradient * @step + 0.5 * @hessian.inverse.vt_self_v(@step,@step)

		# Apply step
		@vector_p = @vector
		@vector = fix_periodicity(@vector + @step)

		# Success
		return :continue
	end
	
end
