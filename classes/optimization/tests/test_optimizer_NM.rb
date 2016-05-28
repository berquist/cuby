#!/usr/bin/env ruby
$:.unshift(File.dirname(__FILE__)+"/../../..").uniq! # Magic line to add path to the executable to $RUBYLIB
require "classes/cuby.rb" # Load Cuby classes and libraries
Cuby.load_framework

require "classes/optimization/optimizer_LBFGS.rb"
require "classes/optimization/optimizer_NM.rb"

#===============================================================================
# Test problem: multidimensional parabola
#===============================================================================

class OptProblemParabola

	def initialize(size)
		@size = size
	end

	def size
		return @size
	end

	def start_point
		a = []
		size.times{ a << 1.0 }
		return Vector.from_array(a)
	end

	def energy(vec)
		sum = 0.0
		vec.each{|x|
			sum += x**2
		}
		return sum
	end

	def grad(vec)
		g = Vector.zero(@size)
		vec.each_with_index{|x,i|
			g[i] = 2.0 * x
		}
		return g
	end
end

#===============================================================================
# Test problem: Rosenbrock function
#===============================================================================

class OptProblemRosenbrock
	A = 1.0
	B = 100.0

	def initialize
		@size = 2
	end

	def size
		return @size
	end

	def start_point
		return Vector[-3,-4]
	end

	def energy(vec)
		return (A - vec[0])**2 + B*(vec[1] - vec[0]**2)**2
	end

	def grad(vec)
		x,y, = vec.to_a
		g = Vector.zero(@size)
		g[0] = -4.0 * B * x * (y-x**2) - 2.0 * (A - x)
		g[1] = 2.0 * B * (y - x**2)
		return g
	end
end

#===============================================================================
# Optimizer test
#===============================================================================

settings = Settings.new
optimizer = OptimizerNM.new(settings)
problem = OptProblemRosenbrock.new
#problem = OptProblemParabola.new(2)

optimizer.starting_vector = problem.start_point

myvector = nil
optimizer.calculation_queue_proc = Proc.new{|calc_id, vector|
	myvector = vector
}

optimizer.calculation_run_proc = Proc.new{|calc_id|
	result = []
	result << problem.energy(myvector)
	if optimizer.what_to_calculate.include?(:gradient)
		result << problem.grad(myvector)
	end
	# Pass results to the optimizer
	result
}

optimizer.print_step_proc = Proc.new{
	puts "--------------------------------------------------------------------------------"
	puts "Cycle: #{optimizer.cycle}"
	puts "Energy: #{optimizer.energy}"
	puts "Point:"
	optimizer.vector.each{|x|
		printf("  %15.6f\n",x)
	}
	$stdout.flush
}

optimizer.check_convergence_proc = Proc.new{|delta_e|
	# Convergence: gradient
	# optimizer.gradient.abs < 1.0e-4

	# Convergence: just energy	
	delta_e.abs < 1.0e-4
}

# Run the optimizer
opt_status = optimizer.run

puts
puts "Status: #{opt_status}"

