
module ProtocolTest
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :tools
	PROTOCOL_DESCRIPTION = "Test calculation against a known result"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation, :required, "Setup for the calculation"]
	]
	#=======================================================================


	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Calculation test')
	end

	def prepare_protocol
		unless @settings.has_block?(:calculation)
			Cuby::error("Input block 'calculation' has to be defined")
		end

		# Create child job
		settings = @settings.block(:calculation)
		# Copy common settings if present but keep what is set in this step
		#!# Copy from settings.root
		#--------------------------------------------------
		# if @settings.has_block?(:calculation_common)
		# 	settings.copy_from!(@settings.block(:calculation_common), :keep)
		# end
		#-------------------------------------------------- 

		@child_job = Job.new(@name, settings)

		# Try to prepare the job
		begin
			@child_job.prepare
		rescue Exception => e
			@failed = true
			@message = e.message
		end
	end

	def run_protocol(queue)
		return if @failed
		begin
			@child_job.run
		rescue Exception => e 
			@failed = true
			@message = e.message
		end
	end

	def cleanup_protocol
		begin
			@child_job.cleanup
		rescue
		end
	end

	def print_protocol
		# Print original results if requested
		if @settings[:test_print_results]
			@child_job.print_protocol
			Cuby::log.puts_v(:normal)
		end

		# Test name to be printed
		if @settings.set?(:test_name)
			testname = sprintf("%-60s", @settings[:test_name])
		else
			testname = ""
		end

		case evaluate_test
		when :ok
			Cuby::log.puts testname + "[   OK    ]"
		when :differs
			lines = @message.split(/\n/)
			Cuby::log.puts testname + "[ DIFFERS ]   " + lines[0]
			(lines.size-1).times{|i|
				Cuby::log.puts "              " + lines[i+1]
			}
		when :failed
			lines = @message.split(/\n/)
			Cuby::log.puts testname + "[ FAILED  ]   " + lines[0]
			(lines.size-1).times{|i|
				Cuby::log.puts "              " + lines[i+1]
			}
		end
		Cuby::log.puts_v(:normal)
	end

	def results
		# Return results hash
		return @child_job.results
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def evaluate_test
		# If calculation failed, return the status, the message already
		# contains the exception message
		if @failed
			return :failed
		end

		case @settings[:test_result_type]
		when :energy
			target = @settings[:test_result].to_f
			difference = @child_job.results.energy - target
			if difference.abs < @settings[:test_threshold]
				return :ok
			else
				@message = "The energy differs by #{difference} kcal/mol"
			       	if target != 0.0
					percent = difference / target.abs * 100
					@message += " (#{percent} %)"
				end
				return :differs
			end
		when :gradient
			target = @settings[:test_result].to_f
			difference = @child_job.results.gradient.rms - target
			if difference.abs < @settings[:test_threshold]
				return :ok
			else
				@message = "RMS of the gradient differs by #{difference} kcal/mol/A"
				return :differs
			end
		when :dataset
			target = @settings[:test_result].to_f
			difference = @child_job.single_number_result - target
			if difference.abs < @settings[:test_threshold]
				return :ok
			else
				@message = "RMSE in the dataset differs by #{difference} kcal/mol"
				return :differs
			end
		end
	end

end

