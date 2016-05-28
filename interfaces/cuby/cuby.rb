module InterfaceCuby
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = ""
	# Interface information
	INTERFACE = :calculation_external
	CAPABILITIES = [:energy, :ghost_atoms]
	MODIFIER = false
	DIRECTORY = "CUBY"
	METHODS = {
		:proxy => [:energy, :ghost_atoms]
	}
	#=======================================================================
	
	CUBY_INP = "cuby_input.yaml"
	CUBY_OUT = "cuby_output.yaml"


	def prepare_interface
		# Prepare calculation directory
		if calc_dir_mkdir(CUBY_INP) == :old
			calc_using_input
			return
		end

		# Make a copy of the settings, copy without overwrite
		calc_settings = @settings.blocks[:calculation]
		# Save geometry in the settings
		calc_settings[:geometry] = @geometry.to_s
		# Copy charge and multiplicity
		calc_settings[:charge] = @settings[:charge]
		calc_settings[:multiplicity] = @settings[:multiplicity]

		# What to calculate
		if @what == [:energy]
			calc_settings[:job] = :energy
		else
			raise "Only energy calculations possible now"
		end

		# Ask for serialized results
		calc_settings[:write_results_yaml] = CUBY_OUT

		# Write settings as input file
		File.open(in_calc_dir(CUBY_INP), "w+") {|f|
			f.puts(calc_settings.to_hash.to_yaml)
		}

		# Save info on the system upon writing the input
		calc_writing_input
	end

	def calculate_interface
		cuby_run_calculation unless @settings[:existing_calc_dir] == :read_results && 
			(FileTest.exist?(in_calc_dir(CUBY_OUT)) || FileTest.exist?(in_calc_dir("RESULTS/"+CUBY_OUT)))
		results = cuby_read_results
		return results
	end

	def cleanup_interface
		calc_dir_delete if @settings[:job_cleanup]
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def cuby_run_calculation
		command = "cd #{calc_dir};"
		command << " cuby4 #{CUBY_INP} 2> cuby.err > cuby.log"
		system(command)
	end

	def cuby_read_results
		if FileTest.exist?(in_calc_dir(CUBY_OUT))
			f = File.open(in_calc_dir(CUBY_OUT),"r")
		elsif FileTest.exist?(in_calc_dir("RESULTS/"+CUBY_OUT))
			f = File.open(in_calc_dir("RESULTS/"+CUBY_OUT),"r")
		else
			Cuby::error "Output from external cuby calculation (#{CUBY_OUT}) not found in\n#{calc_dir}"
		end

		results = YAML::load(f)
		f.close
		return results
	end

end
