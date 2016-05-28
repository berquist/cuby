require "classes/misc/grep.rb"

module ProtocolGaussianDriver
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :driver
	PROTOCOL_DESCRIPTION = "Uses Gaussian as extrenal driver"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation,          :required, "Setup for the calculation called by Gaussian"]
	]
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Cuby calculation with Gaussian driver')
	end

	def prepare_protocol
		# Read geometry
		load_geometry
		# Write input for gaussian
		write_gaussian_driver_input
		# Write input for child cuby calculation
		write_cuby_input
	end

	def run_protocol(queue)
		Cuby::log.puts_v_only(:normal,"Gaussian calculation running...")
		# Execute gaussian driver
		command = ""
		command << "#{@settings[:gaussian_driver_exe]} #{@name}_gaussian_input.com 2> driver.err > driver.out"
		system(command)

		# Check for errors in the output
		if Grep::file_unicode("#{@name}_gaussian_input.log", /Error termination/)
			Cuby::error("Gaussian run failed, see its output for details")
		end
	end

	def print_protocol
		# Print results
		Cuby::log.puts_v_only(:normal,"Gaussian calculation finished")
	end

	def cleanup_protocol
		# Remove unnecessary files
		if @settings[:job_cleanup]
			`rm -f child_cuby_input.yaml`
			`rm -f child_cuby.LOG`
			`rm -f temporary_geo.xyz`
			`rm -f driver.err`
			`rm -f driver.out`
			`rm -f job_gaussian_input.com`
			`rm -f output_for_gaussian.txt`
		end
	end

	def results
		return nil
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def write_gaussian_driver_input
		f = File.open("#{@name}_gaussian_input.com", "w+")

		f.puts "#External=\"#{protocol_dir}/cuby_for_gaussian\" #{@settings[:gaussian_driver_job]}"
		f.puts
		f.puts "Gaussian as a driver for cuby"
		f.puts
		f.puts "0 1" # Charge and multiplicity - does not matter
		# Geometry:
		@geometry.write_xyz(:file => f, :no_header => true, :format => "%25.15f")
		f.puts
		f.puts

		f.close
	end

	def write_cuby_input
		geometry_fn = "temporary_geo.xyz"

		# Update child settings
		settings = @settings.block(:calculation)
		settings[:geometry] = geometry_fn
		settings[:gaussian_external_output] =  "output_for_gaussian.txt"

		# Write the child input file
		settings.save_to_file("child_cuby_input.yaml")
	end


end
