
class Job
	attr_reader :name, :settings
	attr_reader :protocol_name

	#=======================================================================
	# Build list of protocols available in the protocols directory
	#=======================================================================
	@@protocols = {} # List of available protocols
	@@protocols_path = Cuby.install_dir + "/protocols"
	# Fill the list
	Dir.entries(@@protocols_path).sort.each{|dir|
		next if dir == '.' || dir == '..'
		next unless FileTest.directory?(@@protocols_path + "/" + dir)
		@@protocols[dir.downcase.to_sym] = @@protocols_path + "/" + dir
	}

	#=======================================================================
	# Access to the list of protocols
	#=======================================================================
	def Job.available_protocols
		return @@protocols
	end

	#=======================================================================
	# Core functionality
	#=======================================================================
	def initialize(name, settings)
		# Create new job, loading the appropriate protocol module
		@name = name
		@settings = settings

		unless Job.available_protocols.has_key?(@settings[:job])
			Cuby::error("Job type '#{@settings[:job]}' not found.")
		end

		# Load the protocol file and extend this job with the prtotocol module
		protocol_name = @settings[:job].to_s
		protocol_path = Job.available_protocols[@settings[:job]]
		require protocol_path + "/" + protocol_name + ".rb"
		protocol_module = eval("Protocol" + protocol_name.to_s.split("_").map{|s| s.capitalize}.join(''))
		extend protocol_module

		# Save the protocol name
		@protocol_name = @settings[:job]

		# No results yet
		@results = nil
	end

	def print_header
		# Print ASCIIart header for selected calculation
		if @settings.set?(:job_title)
			Cuby.log.print_cuby_logo(:comment => false, :text => @settings[:job_title])

		else
			print_header_protocol
		end
	end

	def prepare
		#: Preparation of calculations etc
		#: Allowed options:
		#* :silent (true/false) - disable printing of the header
		#!# options might be passed through settings
		prepare_protocol
		if @settings[:prepare_only]
			Cuby.log.puts('Calculations have been prepared, calculation will be skipped');
		end
	end

	def run(queue = Cuby.queue)
		return if @settings[:prepare_only]
		# Invalidate old results, we want to calculate it anew
		@results = nil
		# Runs the job
		run_protocol(queue)
	end

	def print
		return if @settings[:prepare_only]
		# Print results to log
		print_protocol
		
		# Save serialized results
		if @settings.set?(:write_results_yaml)
			File.open(@settings[:write_results_yaml], "w+"){|f|
				f.puts(results.to_yaml)
			}
		end
	end

	def cleanup
		return if @settings[:prepare_only]
		cleanup_protocol
	end

	def geometry
		# If there is a single geometry the protocol operates on, return it
		if @geometry
			return @geometry
		else
			Cuby::error("Job #{name} does not yield the geometry it operates on")
		end
	end

	#=======================================================================
	# Optional functionality
	#=======================================================================
	
	# Optionally, the protocols might implement following methods:
	#
	# single_line_result
	# returns a string with the bare result - e.g. energy job would return
	# a string "Energy: 1.23456 kcal/mol" (if applicable)
	#
	# single_number_result
	# returns just one number (if applicable)

	#=======================================================================
	# Helpers
	#=======================================================================
	
	# If more code is added to this section, it should be moved to a separate
	# module (e.g. ProtocolCommon) that will be included or used as ancestor
	# of the protocols that use it 
	
	def load_geometry
		# Loads the geometry from :geometry keyword, applying selection
		# from the :selection keyword (everything by default)
	
		# Load geometry
		@geometry = Geometry.load_from_settings(@settings)

		# Track changes to the geometry, to be printed later
		changes = []

		# Set ghost atoms
		if @settings.set?(:ghost_atoms)
			 @geometry.geometry_from_selection(@settings[:ghost_atoms]).each{|atom|
				 atom.properties[:ghost] = true
			 }
			changes << "some atoms changed to ghosts"
		end

		# Load atomic charges into geometry
		if @settings.set?(:geometry_load_charges)
			charges = AtomicCharges.from_file(@settings[:geometry_load_charges])
			charges.write_to_geometry(@geometry)
			Cuby::log.puts_debug("Atomic charges from file #{@settings[:geometry_load_charges]} loaded into geometry")

		end

		# Apply selection
		unless @settings[:selection] == "auto"
			@geometry = @geometry.geometry_from_selection(@settings[:selection])
			changes << "selection applied"
		end

		# Print the geometry if requested
		if @settings[:print].include?(:input_geometry)
			if changes.size > 0
				s = " (" + changes.join(', ') + ")"
			else
				s = ""
			end
			Cuby::log.puts "Input geometry#{s}:"
			Cuby::log.puts @geometry.to_s
			Cuby::log.puts
		end

		# Charge and multiplicity from .xyz
		if @geometry.info[:xyz_remark] && @settings[:geometry_setup_from_file]
			@settings[:charge] = 0
			@settings[:multiplicity] = 1
			@geometry.info[:xyz_remark].strip.split.each{|entry|
				k,v = entry.split("=")
				if k.downcase == "charge"
					@settings[:charge] = v.to_i
				end
				if k.downcase == "multiplicity"
					@settings[:multiplicity] = v.to_i
				end
			}
		end

		# Charge / multiplicity from other files
		if @geometry.info[:geofile_charge] && @settings[:geometry_setup_from_file]
			@settings[:charge] = 0
			@settings[:multiplicity] = 1
			# Charge, test for integer value (in case it is a sum of atomic charges)
			if @geometry.info[:geofile_charge]
				charge = @geometry.info[:geofile_charge].round
				err = (@geometry.info[:geofile_charge] - charge).abs
				if err > 0.1
					"Charge from geometry file (#{@geometry.info[:geofile_charge]}) is far from integer value"
				end
				@settings[:charge] = charge
			end
			# Multiplicity
			if @geometry.info[:geofile_multiplicity]
				@settings[:multiplicity] = @geometry.info[:geofile_multiplicity]
			end
		end

	end

	def protocol_dir
		return Job.available_protocols[@settings[:job]]
	end

end
