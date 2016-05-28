
module ProtocolMultistep
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :workflow
	PROTOCOL_DESCRIPTION = "Chaining of multiple calculations"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation_X,      :required, "Input for step X (X = step name)"],
		InputBlock[:calculation_common, :optional, "Common settings for all steps"]
	]
	#=======================================================================
	
	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Multistep calculation')
	end

	def prepare_protocol
		# Get list of steps
		@step_names = @settings[:steps].map{|x| x.to_s } # Convert symbols back to strings
		Cuby::log.puts_debug("Job steps:\n")

		# Check subsections
		@step_names.each{|name| 
			if @settings.has_block?(("calculation_"+name).to_sym)
				Cuby::log.puts_debug("   #{name} - setup found")
			else
				# Create non-existent subsections
				@settings.new_block(("calculation_"+name).to_sym)
				Cuby::log.puts_debug("   #{name} - block created")
			end
		}

		# Create step jobs
		@step_jobs = []
		@step_names.each_index{|i|
			name = @step_names[i]
			# Build settings
			settings = @settings.block(("calculation_"+name).to_sym)
			# Copy common settings if present but keep what is set in this step
			if @settings.has_block?(:calculation_common)
				settings.copy_from!(@settings.block(:calculation_common), :keep)
			end

			# Create child job
			job = Job.new(@name+"_step_%03d" % (i+1), settings)
			@step_jobs << job
		}

		# Copy step title decoration to childs
		if @settings.set?(:step_title_decoration)
			@step_jobs.each{|job|
				job.settings[:step_title_decoration] = @settings[:step_title_decoration] unless job.settings.set?(:step_title_decoration)
			}
		end
	end

	def run_protocol(queue)
		queued_jobs = []
		queued_names = []
		@skipped_indexes = []
		@step_names.each_index{|i|
			name = @step_names[i]
			job = @step_jobs[i]

			# Skip
			if job.settings.set?(:skip_step_if_file_found)
				if FileTest.exist?(job.settings[:skip_step_if_file_found])
					# Step header
					print_step_header(job, name)
					# Skip note
					Cuby.log.puts_v(:normal, "Step skipped, file #{job.settings[:skip_step_if_file_found]} exists")
					Cuby.log.puts_v(:normal, "")
					@skipped_indexes << i
					next
				end
			end
			if job.settings.set?(:skip_step_if_file_not_found)
				unless FileTest.exist?(job.settings[:skip_step_if_file_not_found])
					# Step header
					print_step_header(job, name)
					# Skip note
					Cuby.log.puts_v(:normal, "Step skipped, file #{job.settings[:skip_step_if_file_not_found]} not found")
					Cuby.log.puts_v(:normal, "")
					@skipped_indexes << i
					next
				end
			end


			if job.settings[:step_queue]
				Cuby::error("Keyword step_queue is not compatible with geometry: previous_step") if job.settings.set?(:geometry) && job.settings[:geometry].downcase == "previous_step"
				# Queue this step
				job.prepare
				job.run
				queued_jobs << job
				queued_names << name
			else
				# Flush queue
				run_queued(queued_jobs, queued_names)
				# Run immediately
				if @settings[:multistep_print]
					# Step header
					print_step_header(job, name)
				end
				# Update geometry
				if job.settings.set?(:geometry) && job.settings[:geometry].downcase == "previous_step"
					if i == 0
						Cuby::error("Geometry update can not be performed in first step")
					else
						job.settings[:geometry] = @step_jobs[i-1].geometry.to_s
					end

				end
				# Prepare
				job.prepare
				# Run
				job.run
				if @settings[:multistep_print]
					# Job results
					job.print 
				end
			end
		}
		# Run the rest of the queue
		run_queued(queued_jobs, queued_names)
	end

	def cleanup_protocol
		# Cleanup only after everything finished
		@step_jobs.each_index{|i|
			@step_jobs[i].cleanup unless @skipped_indexes.include?(i)
		}
	end

	def print_protocol
		# The final result
		if @settings.set?(:multistep_result_expression)
			Cuby::log.puts_v(:normal, '*' * 80)
			Cuby::log.puts("#{@settings[:multistep_result_name]}: #{single_number_result}")
			Cuby::log.puts_v(:normal, '*' * 80)
			Cuby::log.puts_v(:normal, "")
		end

		# Run a custom code
		if @settings.set?(:multistep_result_eval)
			steps = results_hash
			Cuby::log.puts_v(:normal, '*' * 80)
			Cuby::log.puts_v(:normal, "Processed results:")
			Cuby::log.puts_v(:normal, '*' * 80)
			eval(@settings[:multistep_result_eval])
			Cuby::log.puts_v(:normal, '*' * 80)
			Cuby::log.puts_v(:normal, "")
		end
	end

	def results
		# By default, the result of the last step is returned
		unless @results
			#@results = results_hash[@step_names.last]
			@results = @step_jobs.last.results
		end
		return @results
	end

	#=======================================================================
	# Optional protocol functionality
	#=======================================================================
	
	def single_number_result
		if @settings.set?(:multistep_result_expression)
			steps = results_hash
			value = eval(@settings[:multistep_result_expression])
			Cuby::log.puts_debug "#{@settings[:multistep_result_name]}: #{value.to_s}"
			return value
		else
			Cuby::error("Multistep protocol can not produce single number result\nunless multistep_result_expression keyword is set")
		end
	end

	#=======================================================================
	# Private
	#=======================================================================
	
	def results_hash
		# Optionally, hash of he results can be requested
		unless @results_hash
			@results_hash = {}
			@step_names.each_index{|i|
				@results_hash[@step_names[i]] = @step_jobs[i].results unless @skipped_indexes.include?(i)
			}
		end
		# Return results hash
		return @results_hash
	end

	def print_step_header(job, name)
		decoration_char = ''
		line1 = line2 = false
		if job.settings[:step_title_decoration] != ''
			decoration_char = job.settings[:step_title_decoration][0..0] 
			line1 = line2 = true
			if job.settings[:step_title_decoration].size > 1 && job.settings[:step_title_decoration][1..1] == " "
				line2 = false
			end
		end
		Cuby.log.puts_v(:normal, decoration_char * 80) if line1
		if job.settings.set?(:step_title)
			Cuby.log.puts_v(:normal, job.settings[:step_title])
		else
			Cuby.log.puts_v(:normal, "Step #{name}")
		end
		Cuby.log.puts_v(:normal, decoration_char * 80) if line2
	end

	def run_queued(queued_jobs,  queued_names)
		# Execute the queue
		Cuby.queue.execute
		# Print steps
		queued_jobs.each_index{|j|
			if @settings[:multistep_print]
				# Step header
				print_step_header(queued_jobs[j], queued_names[j])
				# Job results
				queued_jobs[j].print 
			end
		}
		# Empty the list
		queued_jobs.clear
		queued_names.clear
	end
end

