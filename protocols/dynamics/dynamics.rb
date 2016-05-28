require "protocols/dynamics/molecular_dynamics.rb"
require "protocols/dynamics/trajectory.rb"

# TBD:
# keywords for trajectory writing
# velocities: reading
# solve sharing keywords with optimization
# test everything

module ProtocolDynamics
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :driver
	PROTOCOL_DESCRIPTION = "Molecular dynamics engine"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:trajectory_2, :optional, "Setup for the optional secondary trajectory file"]
	]
	#=======================================================================
	
	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Molecular dynamics')
	end

	def prepare_protocol
		# Read geometry
		load_geometry

		# Create calculation
		@calculation = Calculation.new(@name, @settings, @geometry)
		@calculation.parent_protocol = self

		# Prepare calculation
		@calculation.prepare([:energy, :gradient])

		# Create and setup the MD object
		@molecular_dynamics = MolecularDynamics.new(@geometry, @settings)

	end

	def run_protocol(queue)
		# Initialize trajectory writer
		trajectory_writers = [Trajectory::Writer.new(@settings, @molecular_dynamics)]

		if @settings.has_block?(:trajectory_2)
			@settings.block(:trajectory_2).input_file = @settings.input_file
			trajectory_writers << Trajectory::Writer.new(@settings.block(:trajectory_2), @molecular_dynamics, 2)
		end



		# Printing proc
		@time_last = Time.now
		@molecular_dynamics.print_step_proc = Proc.new{
			# Timing
			newtime = Time.now
			Cuby::log.puts_debug("MD step took #{'%.2f' % (newtime - @time_last)} seconds")
			@time_last = newtime

			puts "--------------------------------------------------------------------------------"
			puts "Cycle: #{@molecular_dynamics.cycle}"
			puts "Energy: #{@molecular_dynamics.energy}"
			puts "Kinetic energy: #{@molecular_dynamics.kinetic_energy}"
			puts "Temperature: #{@molecular_dynamics.temperature}"
			if @settings.set?(:temperature_target)
				puts "Thermostat: #{@molecular_dynamics.thermostat_temp}"
			end
			# Print energy components
			@calculation.results.print_energy_decomposition if @settings[:print].include?(:energy_decomposition)
			@calculation.results.print_dipole if @settings[:print].include?(:dipole)
			@calculation.results.print_gradient if @settings[:print].include?(:gradient)
			# Write trajectory
			trajectory_writers.each{|writer| writer.write}
		}


		# Run simulation
		@molecular_dynamics.run{
			@calculation.send_to_queue(queue)

			# Pass results to the optimizer
			[
				@calculation.results.energy,
				@calculation.results.gradient
			]

		}

		# Print final status
		puts "================================================================================"
		puts "Simulation finished"

		# Close the trajectory
		trajectory_writers.each{|writer| writer.write_last}

		# Write the last geometry as a restart file
		@geometry.write_xyz(:file => "last.xyz", :velocities => @molecular_dynamics.velocities_full)
	end

	def print_protocol
	end

	def cleanup_protocol
		@calculation.cleanup
	end

	#-----------------------------------------------------------------------
	# Expose the MD data
	#-----------------------------------------------------------------------
	
	def molecular_dynamics
		return @molecular_dynamics
	end

end
