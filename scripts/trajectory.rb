#!/usr/bin/env ruby
$:.unshift(File.dirname(__FILE__)+"/..").uniq! # Magic line to add path to the executable to $RUBYLIB
require "classes/cuby.rb" # Load Cuby classes and libraries

class TrajectoryAction
	def initialize(settings)
		@counter = 0
	end

	def step(geometry)
		@counter += 1
	end

	def finish
	end
end

class TrajectoryActionSort < TrajectoryAction

	def initialize(settings)
		super(settings)
		@frames = []
	end

	def step(geometry)
		if geometry.info[:xyz_remark]
			@frames << geometry.deep_copy
		else
			$stderr.puts "Remark line not found in frame #{@counter}, frame skipped"
		end
		super(geometry)
	end

	def finish
		# Get energies
		@frames.each{|geo|
			geo.info[:energy] = geo.info[:xyz_remark].split[1].to_f
		}
		# Sort
		@frames.sort!{|a,b| a.info[:energy] <=> b.info[:energy]}
		#Write sorted
		@frames.each{|geo|
			geo.write_xyz(:second_line => geo.info[:xyz_remark])
		}
	end
end

class TrajectoryActionRmsdFilter < TrajectoryAction

	def initialize(settings)
		super(settings)
		@frames = []
		@cutoff = 1.5
	end

	def step(geometry)
		# Check RMSD against existing frames
		$stderr.puts "#{@counter}" if @counter % 10 == 0
		skip = false
		@frames.each{|frame|
			fitted = geometry.fit_rmsd(frame)
			rmsd = fitted.rmsd(frame)
			if rmsd < @cutoff
				skip = true
				break
			end
			# Try mirror image
			geometry.each{|atom| atom.x = -atom.x}
			fitted = geometry.fit_rmsd(frame)
			rmsd = fitted.rmsd(frame)
			if rmsd < @cutoff
				skip = true
				break
			end
		}
		@frames << geometry.deep_copy unless skip
		super(geometry)
	end

	def finish
		# Write selected frames
		@frames.each{|geo|
			geo.write_xyz(:second_line => geo.info[:xyz_remark])
		}
	end
end


class TrajectoryActionAverage < TrajectoryAction

	def step(geometry)
		if @counter == 0
			@sum = geometry.to_vector
			@ref_geo = geometry.deep_copy
		else
			@sum += geometry.to_vector
		end

		super(geometry)
	end

	def finish
		@ref_geo.update_from_vector!(@sum / @counter)
		@ref_geo.write_xyz
	end
end

class TrajectoryActionFrameDirs < TrajectoryAction

	def initialize(settings, limit)
		super(settings)
		@limit = limit
	end

	def step(geometry)
		if @counter < @limit
			dirname = sprintf("frame_%05d",@counter)
			system "mkdir #{dirname}"
			geometry.write_xyz(:file => "#{dirname}/geo.xyz", :second_line => geometry.info[:xyz_remark])
		end
		super(geometry)
	end
end

Cuby.application {
	# Load Cuby classes and libraries
	Cuby.load_framework

	# Parse commandline and build settings
	settings = Settings.from_commandline(ARGV, false)

	# Configure Cuby module using the settings
	Cuby.configure(settings)

	# Open the geometry file
	filename = settings[:geometry]
	f = File.open(filename, "r")

	# Create the action
	#action = TrajectoryActionAverage.new(settings)
	#action = TrajectoryActionSort.new(settings)
	#action = TrajectoryActionRmsdFilter.new(settings)
	action = TrajectoryActionFrameDirs.new(settings, 20)

	# Iterate over each frame - for now, only .xyz format
	geo = Geometry.new
	while GeoFile::XYZ.read(geo, {:file => f})
		action.step(geo)
	end

	action.finish
	f.close

	Cuby.finalize # Cleanup etc
}
