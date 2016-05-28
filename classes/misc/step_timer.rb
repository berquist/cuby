# Simple tool for printing the timing of individual steps of the calculation

class StepTimer
	def initialize(settings, verbosity_level = :normal)
		@print = settings[:print].include?(:step_timing)
		@verbosity = verbosity_level
	end

	def step(step_name)
		Cuby::log.puts_v(@verbosity, "Timer: entering #{step_name}") if @print
		start_t = Time.now
		retval = yield
		stop_t = Time.now
		Cuby::log.puts_v(@verbosity, "Timer: #{stop_t - start_t} seconds (step #{step_name})") if @print

		return retval
	end
end

