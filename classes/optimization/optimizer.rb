
class Optimizer
	# Externally provided procs
	attr_accessor :calculation_queue_proc # One argument passed: the current vector
	attr_accessor :calculation_run_proc # Return value: [energy, gradient]
	attr_accessor :check_cart_proc # Return value: new step
	attr_accessor :print_step_proc # No arguments passed
	attr_accessor :check_convergence_proc # One argument: delta_e, Return value: true if converged
	attr_accessor :hessian_estimate_proc # Builds a hessian estimate

	# Readers for optimizer variables
	attr_reader :cycle
	attr_reader :vector
	attr_reader :energy
	attr_reader :gradient
	attr_reader :delta_e

	# Periodicity handling
	# Array of 2 Vectors -> min/max bounds for each coordinate
	# nil -> no periodicity treatment
	attr_accessor :periodicity
	
	# Statistic counters
	attr_accessor :stat_counters

	def initialize (settings)
		@cycle = 0
		@stat_counters = {}
		@stat_counters.default = 0
		@settings = settings
		
		@delta_e = nil

		@periodicity = nil
		
		# Prepare linesearch
		case @settings[:linesearch]
		when :none
			require "classes/optimization/linesearch_none.rb"
			extend LinesearchNone
		when :quadratic_ls
			require "classes/optimization/linesearch_quadratic_ls.rb"
			extend LinesearchQuadraticLS
		when :old
			require "classes/optimization/linesearch_old.rb"
			extend LinesearchOld
		when :cubic
			require "classes/optimization/linesearch_cubic.rb"
			extend LinesearchCubic
		end
		
	end

	def parallel
		# By default, optimizers are not parallel
		return false
	end

	def what_to_calculate
		# By default, optimizer needs energy and gradient. This method can be overriden in
		# ancestor optimizer.
		return [:energy, :gradient]
	end

	def check
		# Check necessary setup before optimize is run
		raise "Calculation_queue proc must be supplied to the Optimizer" unless @calculation_queue_proc
		raise "Calculation_run proc must be supplied to the Optimizer" unless @calculation_run_proc
	end

	def starting_vector=(vec)
		# Setter for @vector, to be used only upon initialization of the optimizer
		@vector = vec
	end
	
	#~ def step
		#~ #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#~ # This is only a template for actual optimizers #!# might be old
		#~ #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#~ raise "Template, should be overriden by actual optimizer"
		#~ 
		#~ # Save current data
		#~ @energy_p = @energy
		#~ @grad_p = @gradient
		#~ 
		#~ # Call the calculation, printing and convergence check
		#~ @energy, @gradient = @calculation_proc.call(@vector)
		#~ 
		#~ # Energy change
		#~ @delta_e = @energy - @energy_p if @energy_p
		#~ 
		#~ # Print step
		#~ @print_step_proc.call
		#~ 
		#~ # Check convergence
		#~ if @cycle > 0
			#~ return :converged if @check_convergence_proc.call(@delta_e) 
		#~ end
		#~ 
		#~ @step = some_step
		#~ @vector += @step
		#~ 
		#~ # Success
		#~ return :continue
	#~ end
	
	def hessian_estimate_initialization
		if hessian_estimate_proc
			Cuby::log.puts_debug("Hessian initialization: calling supplied Proc")
			hessian_estimate_proc.call
		else
			Cuby::log.puts_debug("Hessian initialization: calling default method")
			hessian_estimate_default
		end
	end

	def hessian_estimate_default
		# Do nothing - default for optimizers that do not need hessian initialization
	end
	
	def run
		# Call hessian initialization if available
		hessian_estimate_initialization

		step_status = :continue
		while @cycle < @settings[:maxcycles] && step_status == :continue
			step_status = step
			@cycle += 1 if step_status == :continue
		end

		if step_status == :continue && @cycle == @settings[:maxcycles]
			step_status = :maxcycles_exceeded
			@cycle -= 1 # Revert the counter to the last valid cycle
		end

		return step_status
	end

	def print_stats
		return if @stat_counters.size == 0
		Cuby::log.puts "Statistics:"
		@stat_counters.each_pair{|key, value|
			Cuby::log.puts "   #{key}: #{value}"
		}
	end

	#=======================================================================
	# Shared methods
	#=======================================================================
	
	def fix_periodicity(vector)
		return vector unless @periodicity

		newvec = vector.deep_copy
		p_min = @periodicity[0]
		p_max = @periodicity[1]

		count = 0
		vector.each_index{|i|
			period = p_max[i] - p_min[i]
			next if period == 0.0
			newvec[i] -= period while newvec[i] > p_max[i]
			newvec[i] += period while newvec[i] <= p_min[i]
			if newvec[i] != vector[i]
				count += 1
			end
		}
		Cuby::log.puts_debug("Periodicity fix: #{count}")

		return newvec
	end
end
