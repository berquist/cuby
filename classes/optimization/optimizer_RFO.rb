require "classes/optimization/optimizer.rb"
require "classes/optimization/optimizer_common.rb"

class OptimizerRFO < Optimizer
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

	def step_vector_rfo
		# Gradient as matrix
		gmat = @gradient.to_matrix

		# Build extended hessian
		hext = Matrix.zero(@hessian.n + 1)
		hext.paste!(0,0,@hessian)
		hext.paste!(0,@hessian.n,gmat)
		hext.paste!(@hessian.n,0,gmat.transpose)

		# Diagonalize extended hessian
		evec, real, imag = hext.eigensystem

		# Get minimum component and its index
		min, min_i, min_j = real.min_with_index

		a = evec.column_as_array(min_i)
		last = a.pop
		s = Vector.from_array(a) * (1.0/last) 

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
		hessian_update_qn unless step_status == :refused

		# Propose step
		direction = fix_periodicity(step_vector_rfo)
		@step = fix_periodicity(linesearch(direction))

		# Predicted energy change
		@predicted_de = @gradient * @step + 0.5 * @hessian.inverse.vt_self_v(@step,@step)

		# Apply step
		@vector_p = @vector
		@vector = fix_periodicity(@vector + @step)

		# Success
		return :continue
	end

end
