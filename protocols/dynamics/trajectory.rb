module Trajectory
	class Writer
		def initialize(settings, md_object, number = nil)
			# Keep reference to the job
			@md_object = md_object

			# Writing frequency
			@freq = settings[:trajectory_freq]

			# Trajectoty file format
			@format = settings[:trajectory_format] 

			# Are velocities needed?
			if @format == :xyz_velo
				@write_velocities = true
			end

			# Selection saved in trajectory
			case settings[:trajectory_selection]
			when "%all()"
				@selection = :all
			when "md_region"
				if settings[:md_region] == "%all()"
					@selection = :all
				else
					@selection = :region
				end
			else
				@selection = @md_object.geometry_full.atomlist_from_selection(settings[:trajectory_selection])
			end

			# Build filename
			case @format
			when :xyz
				fn_suffix = ".xyz"
			when :xyz_velo
				fn_suffix = ".xyz"
			end

			# Secondary trajectory option
			if number
				fn_suffix = "_#{number}" + fn_suffix
			end

			if settings.set?(:trajectory_file)
				filename = settings[:trajectory_file]
			else
				if settings.input_file
					filename = "trajectory_" + File.basename(settings.input_file).gsub(/\.[^\.]*$/,'') + fn_suffix
				else
					filename = "trajectory" + fn_suffix
				end
			end

			# Open the file
			@file = File.open(filename, "w+")

		end

		def writestep
			# Header line, if applicable
			header = "Cycle: #{@md_object.cycle}"
			header << " E: #{@md_object.energy} kcal/mol"
			header << " E_kin: #{@md_object.kinetic_energy} kcal/mol"
			header << " T: #{@md_object.temperature} K"

			# Geometry and associated velocities
			case @selection
			when :all
				geometry = @md_object.geometry
				velocities = @md_object.velocities_full if @write_velocities
			when :region
				geometry = @md_object.geometry
				velocities = @md_object.velocities if @write_velocities
			else
				geometry = @md_object.geometry_full.geometry_from_list(@selection)
				if @write_velocities
					velo = @md_object.velocities_full if @write_velocities
					velocities = []
					@selection.each{|i| velocities << velo[i]}
				end
			end

			# Write
			case @format
			when :xyz
				geometry.write_xyz(:file => @file, :second_line => header)
			when :xyz_velo
				geometry.write_xyz(:file => @file, :second_line => header, :velocities => velocities)
			end
			@file.flush
		end

		def write
			return if @freq == 0
			# Write each @freq cycles
			if @md_object.cycle % @freq == 0
				writestep
				@last_cycle_written = @md_object.cycle
			end
		end

		def write_last
			return if @freq == 0
			# Write final geometry if it was not written yet
			unless @last_cycle_written == @md_object.cycle
				writestep
			end
		end

		def close
			@file.close
		end
	end
end
