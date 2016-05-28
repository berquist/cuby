require "classes/geometry/geometry_series.rb"

module ProtocolScan
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :batch
	PROTOCOL_DESCRIPTION = "Scan over a series of geometries"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation, :required, "Setup for the calculation"]
	]
	#=======================================================================


	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Scan')
	end

	def prepare_protocol
		# Generate geometries
		unless @settings[:scan_generator] == :none
			# Load the initial geometry
			load_geometry

			# Initialize the generator
			case @settings[:scan_generator]
			when :grid_xyz
				require "protocols/scan/classes/scan_generator_grid_xyz"
				generator = ScanGeneratorGridXyz.new(@settings)
			end

			# Run the generator
			if @settings[:scan_generator_batchsize] == 0
				generator.run(@geometry, "generated_scan.xyz")
			else
				dirlist = generator.build_batches(@geometry, @settings[:scan_generator_batchsize], "scan_batch", "generated_scan.xyz")
				# Modify settings
				@settings[:geometry] = "generated_scan.xyz"
				@settings[:scan_generator] = :none
				# Save the new settings in each generated directory
				dirlist.each{|dirname|
					File.open(dirname + "/input.yaml", "w+"){|f|
						f.puts @settings.to_hash.to_yaml
					}
				}
				Cuby::log.puts "Batch calculations prepared 'scan_batch_...' directories, exitting"
				exit
			end

			@settings[:geometry] = "generated_scan.xyz"

		end

		# Check whether geometry is an .xyz file
		unless FileTest.file?(File.expand_path(@settings[:geometry]))
			Cuby::error("For the scan job, geometry keyword must reference to a .xyz file")
		end

		# Check the presence of the Calculation block
		unless @settings.has_block?(:calculation)
			Cuby::error("Block 'calculation' not found in the input")
		end

		# Determine the mode of operation and load the module implementing it
		case @settings[:scan_mode]
		when :serial
			self.extend ProtocolScanSerial
		when :parallel
			self.extend ProtocolScanFullParallel
		end

		# Run the mode-specific prepare
		prepare_protocol2
	end

end

module ProtocolScanSerial

	def prepare_protocol2
		# Load the first geometry or a template
		g = Geometry.new
		if @settings.set?(:geometry_template)
			g.read_file(File.expand_path(@settings[:geometry_template]))
		else
			f = File.open(File.expand_path(@settings[:geometry]))
			g.read_xyz(:file => f)
			f.close
		end


		# Pass it to the calculation block
		calc_settings = @settings.block(:calculation)
		calc_settings[:geometry] = g.to_s
		
		# Set up and prepare the child calculation
		@child_job = Job.new(@name + "_calc", calc_settings)
		Cuby::error("The serial scan is compatible only with protocols that support 'update_geometry' method") unless @child_job.respond_to?(:update_geometry)
		@child_job.prepare

	end

	def run_protocol(queue)
		# Empty results array
		@results = []
		# Iterate over the geometries
		f = File.open(File.expand_path(@settings[:geometry]))
		g = Geometry.new
		while g.read_xyz(:file => f)
			# Get label from second line of the xyz file
			label = g.info[:xyz_remark]
			if label =~ /^scan:/
				label.gsub!(/^scan:\s*/,"")
			else
				label = nil
			end

			@child_job.update_geometry(g)

			# Run the calculation
			success = true
			begin
				@child_job.run
				# Store results
				@results << @child_job.results
			rescue
				success = false
				# Rebuild the calculation
				@child_job.cleanup
				@child_job.prepare
			end

		
			# Print	
			if @settings[:scan_print_mode] == :numbers && @child_job.respond_to?(:single_number_result)
				# simple oneline printing
				if label 
					s = label + "\t"
				else
					s = ""
				end
				if success
					s += "%20.6f" % @child_job.single_number_result
				else
					s += "%20s" % "failed"
				end
				Cuby::log.puts(s)
			elsif @settings[:scan_print_mode] == :simple && @child_job.respond_to?(:single_line_result)
				# simple oneline printing
				if label 
					s = label + "\t"
				else
					s = ""
				end
				if success
					s += "%20s" % @child_job.single_line_result.to_s
				else
					s += "%20s" % "failed"
				end
				Cuby::log.puts(s)
			else
				Cuby::log.puts("Geometry label: #{label}") if label
				if success
					@child_job.print_protocol
				else
					Cuby::log.puts("Calculation of this point failed")
				end
			end

			# Clean the geometry
			g = Geometry.new
		end
		f.close
	end

	def cleanup_protocol
		# Clean child job
		@child_job.cleanup
	end

	def print_protocol
	end

	def results
		return @results
	end
end

module ProtocolScanFullParallel

	def prepare_protocol2
		geometries = GeometrySeries.new(settings, "scan:")
		@item_labels = geometries.labels

		Cuby::error("Scan: no geometries read") if geometries.size == 0
		Cuby::warning("Scan: only one geometry read") if geometries.size == 1

		Cuby.log.puts_debug("Read geometries.size geometries")

		unless @settings.has_block?(:calculation)
			Cuby::error("Input for a scan job must contain subsection named 'calculation' which contains the calculation setup")
		end
		unless @settings.set?(:calculation,:job)
			Cuby::error("Job type must be defined in the 'calculation' subsection of the input")
		end

		# Create child jobs
		@item_jobs = []
		geometries.each_index{|i|
			name = sprintf("frame%03d",i)
			blockname = ("calculation_" + name).to_sym
			@settings.new_block(blockname)
			@settings.block(blockname).copy_from!(@settings.block(:calculation), :keep, [/^calculation_frame/])
			settings = @settings.block(blockname)
			settings[:geometry] = geometries[i].to_s
			@item_jobs[i] = Job.new(@name + "_" + name, settings)
		}

		# Prepare child jobs
		@item_jobs.each{|job| job.prepare }

	end

	def run_protocol(queue)
		# Submit jobs to queue
		@item_jobs.each{|job| job.run}
	end

	def cleanup_protocol
		# Clean child jobs
		@item_jobs.each{|job| job.cleanup}
	end

	def print_protocol
		# Print results
		if @settings[:scan_print_mode] == :numbers && @item_jobs.first.respond_to?(:single_number_result)
			# simple oneline printing
			@item_jobs.each_index{|i|
				if @item_labels.size > 0
					s = @item_labels[i] + "\t" 
				else
					s = ""
				end
				s += "%20.6f" % @item_jobs[i].single_number_result
				Cuby::log.puts(s)
			}
		elsif @settings[:scan_print_mode] == :simple && @item_jobs.first.respond_to?(:single_line_result)
			# simple oneline printing
			@item_jobs.each_index{|i|
				if @item_labels.size > 0
					s = @item_labels[i] + "\t" 
				else
					s = ""
				end
				s += @item_jobs[i].single_line_result
				Cuby::log.puts(s)
			}
		else
			# printing using child protocol's print
			@item_jobs.each_index{|i|
				Cuby::log.puts("Geometry label: #{@item_labels[i]}") if @item_labels.size > 0
				@item_jobs[i].print_protocol
			}
		end
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			# Get results
			@results = []
			@item_jobs.each{|job|
				@results << job.results
			}
		end
		return @results
	end
end
