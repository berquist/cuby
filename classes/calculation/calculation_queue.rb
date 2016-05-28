require 'thread'

class CalculationQueue
	def initialize(threads)
		# Initialize empty queue
		@calculations = []

		# Load serial or parallel queue as requested
		if threads == 1
			Cuby::log.puts_debug "Serial queue created"
			extend CalculationQueueSerial
		else
			Cuby::log.puts_debug "Parallel queue created using #{threads} threads"
			extend CalculationQueueParallel
			@threads = threads
		end
	end

	def add_to_queue(calculation)
		Cuby::log.puts_debug "Calculation \"#{calculation.name}\" added to queue"
		@calculations << calculation
	end

	def sort_largest_first!
		@calculations.sort!{|a, b|
			a.priority <=> b.priority
		}
		@calculations.reverse!
		return nil
	end
end

module CalculationQueueSerial
	def execute
		if @calculations.size == 0
			Cuby::log.puts_debug "Queue executed but is empty"
			return 0 
		end

		sort_largest_first!

		size = @calculations.size
		Cuby::log.puts_debug "Queue running - #{size} calculations"
		while calc = @calculations.shift
			Cuby::log.puts_debug "Calculating \"#{calc.name}\" (priority #{calc.priority})"
			calc.calculate(0)
		end
		Cuby::log.puts_debug "Queue execution finished"
		return size
	end
end

module CalculationQueueParallel
	def execute
		if @calculations.size == 0
			Cuby::log.puts_debug "Queue executed but is empty"
			return 0 
		end

		sort_largest_first!

		size = @calculations.size
		Cuby::log.puts_debug "Queue running - #{size} calculations, running #{@threads} in parallel "
		semaphore = Mutex.new
		index = 0
		controllers = []
		@threads.times {|i|
			controllers[i] = Thread.new {
				semaphore.lock
				while index < @calculations.size
					this_i = index
					index += 1
					Cuby::log.puts_debug "   Calculating \"#{@calculations[this_i].name}\" (priority #{@calculations[this_i].priority})"
					semaphore.unlock
					@calculations[this_i].calculate(i)
					semaphore.lock
				end
				semaphore.unlock
			}
			controllers[i].abort_on_exception = true
		}
		controllers.each_index {|i| controllers[i].join}

		@calculations = []
		Cuby::log.puts_debug "Queue execution finished"

		return size
	end
end
