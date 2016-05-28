require "erb"

class QueueSystem
	def initialize(settings, arguments)
		@settings = settings

		# Hard-coded defaults
		@submit_script = "_queue_submit.sh"
		@run_script = "_queue_run.sh"
		@shell = "/bin/bash"

		# Check if job is already submitted
		if FileTest.exist?(@submit_script)
			$stderr.puts "Job was already submitted."
			$stderr.puts "To submit it again, delete the _queue* scripts first."
			$stderr.puts
			exit
		end

		# Select queue system module
		case @settings[:queue_system]
		when :run_in_background
			extend QueueSystemBackground
		when :pbs
			extend QueueSystemPbs
		when :sge
			extend QueueSystemSge
		end

		# Process arguments - remove queuing option
		@arguments = []
		arguments = arguments.dup
		while a = arguments.shift
			if a == "-q" || a.downcase == '--queue-submit' || a.downcase == '--queue_submit'
				# Remove also the value
				arguments.shift
			else
				@arguments << a
			end
		end
	end

	def submit!
		# Prepare the scripts
		prepare_run_script
		prepare_submit_script

		# Make them executable
		system("chmod u+x #{@run_script}")
		system("chmod u+x #{@submit_script}")

		# Run the submit script
		run_submit_script
	end

	#=======================================================================
	# Job parameter accessors
	#=======================================================================
	
	def queue_name
		# Name of the queue used for submission
		# This allows to set it either from input or from config file
		if @settings[:queue_name] == "N/A"
			Cuby::error("Keyword 'queue_name' not set")
		else
			return @settings[:queue_name]
		end
	end

	def job_name
		if @settings.set?(:queue_jobname)
			s = @settings[:queue_jobname]
		elsif @settings.set?(:job_title)
			# Use job title if queue_jobname is not present
			s = @settings[:job_title]
		else
			s = @settings[:interface].to_s
		end
		# Remove spaces
		s.gsub!(/\s+/,"_")
		# Should not start with number - prepend j if needed
		s = "j" + s if s =~ /^[0-9]/

		return s
	end

	def no_of_cpus
		# Determines the number of CPUs to be alocated
		if @settings.set?(:queue_parallel)
			# If queue keyword is set, use that
			return @settings[:queue_parallel]
		else
			# Otherwise get it from calculation setup
			return @settings[:cuby_threads] * @settings[:parallel]
		end
	end

	def scratch_dir
		# Scratch directory, witch check
		# This allows to set it either from input or from config file
		if @settings[:queue_scratch_dir] == "N/A"
			Cuby::error("Keyword 'queue_scratch_dir' not set")
		else
			return  @settings[:queue_scratch_dir]
		end
	end

	def sge_pe
		# SGE only - parallel environment option
		# use default - shared memory
		return "shm"
	end
end

#===============================================================================
# Sample module: job executed in background
#===============================================================================

module QueueSystemBackground
	# This is an example of simple job submission - the job is run in background

	def prepare_run_script
		f = File.open(@run_script, "w+")
		f.puts "#!#{@shell}"
		f.puts "cuby4 -q no #{@arguments.join(' ')} &> #{@settings[:queue_output_file]}"
		f.close
	end

	def prepare_submit_script
		f = File.open(@submit_script, "w+")
		f.puts "#!#{@shell}"
		f.puts "./#{@run_script} &"
		f.close
	end

	def run_submit_script
		# Message
		puts "Cuby job has been sent to background"
		# Run the submit script
		system("./#{@submit_script}")
	end

end

module QueueSystemPbs

	def prepare_run_script
		# Write script
		f = File.open(@run_script, "w+")
		# PBE-specific header template
		template_pbe = ERB.new(IO.read(Cuby.install_dir + "/classes/queue_system/templates/run_pbe.erb"))
		f.print(template_pbe.result(binding))
		# Common job handling template
		template_common = ERB.new(IO.read(Cuby.install_dir + "/classes/queue_system/templates/run_common.erb"))
		f.print(template_common.result(binding))
		f.close
	end

	def prepare_submit_script
		# Write file
		f = File.open(@submit_script, "w+")
		f.puts "#!#{@shell}"
		f.puts "#{@settings[:queue_qsub_command]} #{@settings[:queue_qsub_options]} ./#{@run_script}"
		f.close
	end

	def run_submit_script
		# Message
		puts "Cuby job has been submitted via PBS"
		# Run the submit script
		system("./#{@submit_script}")
	end
end

module QueueSystemSge

	def prepare_run_script
		# Write script
		f = File.open(@run_script, "w+")
		# SGE-specific header template
		template_sge = ERB.new(IO.read(Cuby.install_dir + "/classes/queue_system/templates/run_sge.erb"))
		f.print(template_sge.result(binding))
		# Common job handling template
		template_common = ERB.new(IO.read(Cuby.install_dir + "/classes/queue_system/templates/run_common.erb"))
		f.print(template_common.result(binding))
		f.close
	end

	def prepare_submit_script
		# Write file
		f = File.open(@submit_script, "w+")
		f.puts "#!#{@shell}"
		f.puts "#{@settings[:queue_qsub_command]} -cwd -q #{queue_name} -pe #{sge_pe} #{no_of_cpus} -N #{job_name} #{@settings[:queue_qsub_options]} ./#{@run_script}"
		f.close
	end

	def run_submit_script
		# Message
		puts "Cuby job has been submitted via SGE"
		# Run the submit script
		system("./#{@submit_script}")
	end
end
