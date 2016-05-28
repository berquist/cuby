require "classes/optimization/optimizer.rb"
require "classes/optimization/optimizer_common.rb"

class OptimizerLBFGS < Optimizer
	include OptimizerCommon

	# Hessian estimate handling
	ALLOWED_H_ESTIMATE = [:number_inv, :vector]
	
	def initialize(settings)
		# Generic initialization
		super(settings)

		# Trust radius initialization
		@maxstep = @settings[:trust_radius]
		@maxstep_max = @settings[:trust_radius_max]
		@maxstep_min = @settings[:trust_radius_min]

		# l-BFGS setup
		@lbfgs_n = @settings[:lbfgs_n]

		@oldgrad = []
		@oldvec = []
		
		@step_scaling = @settings[:opt_tr_mode]
		@step_scaling = :abs if @step_scaling == :default
	end

	def hessian_estimate_default
		# Default hessian initialization
		@h0 = 1.0
	end
	
	def initialize_hessian(argument)
		# Initialize hessian
		if argument.kind_of?(Numeric)
			# Size provided, use an universal diagonal guess
			@h0 = 1.0 / argument
		elsif argument.kind_of?(Vector)
			# Vector provided, act like diagonal matrix
			@h0 = argument
			# Invert the numbers
			@h0.each_index{|i|
				@h0[i] = 1.0 / @h0[i]
			}
		else
			raise "Argument of initialize_hessian for LBFGS must be either a number or a vector."
		end
	end
	
	def direction(oldgrad, oldvec)
		# record size
		m = oldgrad.size - 1
		m = @lbfgs_n if m > @lbfgs_n
		
		# LBFGS step calculation
		q = oldgrad[0]
		s = []
		y = []
		ro = []
		alpha = []
		i = 1
		while i < m # Iterate over all saved records
			s[i] = fix_periodicity(oldvec[i-1] - oldvec[i])
			y[i] = oldgrad[i-1] - oldgrad[i]
			ro[i] = 1.0 / (y[i].dot(s[i]))
			alpha[i] = ro[i] * (s[i].dot(q))
			
			q -= alpha[i] * y[i]
			
			i += 1
		end
		
		# Update of numerical hessian estimate according to Nocedal 1989
		# Diagonal internal estimate does not need to be updated
		if @cycle >= 2 and @h0.kind_of?(Numeric)
			@h0 = y[1].dot(s[1]) / y[1].abs**2
		elsif @cycle >= 2 and @h0.kind_of?(Vector)
			h0_nocedal = y[1].dot(s[1]) / y[1].abs**2
			# preserves ratio in @h0, torsion's value = Nocedal's estimate
			#!# might be unstable if extreme value arises from some hessian estimates, rather average of few highest values?
			#@h0 *= h0_nocedal / @h0.max
			#!# for computed hessian estimates, a real update would be better
		end
		
		if @h0.kind_of?(Numeric)
			z = q * @h0
		else
			# act like diagonal matrix
			z = @h0.elementwise_multiply(q)
		end
		
		i = m - 1
		while i > 0
			beta = ro[i]* y[i].dot(z)
			z += s[i] * (alpha[i] - beta)
			
			i -= 1
		end
		
		z *= -1.0
		
		return z
	end
	
	def step
		Cuby::log.puts_debug "=" * 80
		
		# Save current data
		@energy_p = @energy
		
		# Call the calculation, printing and convergence check
		if @reuse_point
			@energy, @gradient = @reuse_point
			@reuse_point = nil
		else
			@real_vector = @calculation_queue_proc.call(0,@vector)
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
		
		# Save coordinates
		@oldgrad.unshift(@gradient)
		@oldvec.unshift(@real_vector)
		#~ @oldvec.unshift(@vector)
		@oldgrad.slice!(@lbfgs_n + 1) # Keep only last N records + 1 so that step can be discarded
		@oldvec.slice!(@lbfgs_n + 1)


		
		# Propose step, update vector
		@step = fix_periodicity(direction(@oldgrad, @oldvec))

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
		# @step = scale_step(direction, @maxstep)
		
		Cuby::log.puts_debug "Step internal: #{@step.max_abs}"
		
		@vector = fix_periodicity(@vector + @step)
		
		# Success
		return :continue
	end
	
end
