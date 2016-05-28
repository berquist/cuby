
# 1) scan the separating ratios
# 2) interpolate optimal ratio
# 2) run constrained and unconstrained calculations

module ProtocolChargeTransfer
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :composite_calculation
	PROTOCOL_DESCRIPTION = "Charge transfer energy in non-covalent dimers"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation, :required, "Setup for the cDFT calculation"]
	]
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Constrained DFT charge transfer')
	end

	def prepare_protocol
		# Input check
		unless @settings.has_block?(:calculation)
			Cuby::error("Input for a scan job must contain subsection named 'calculation' which contains the calculation setup")
		end

		# Read geometry
		load_geometry

		# Population analysis for different positions of the plane
		@scan_ratios = []
		r = @settings[:ct_scan_interval][0]
		while r <= @settings[:ct_scan_interval][2]
			@scan_ratios << r
			r += @settings[:ct_scan_interval][1]
		end

		# Create the calculations
		@calculations = []
		@scan_ratios.each_index{|i|
			name = sprintf("scan%03d",i)
			blockname = ("calculation_" + name).to_sym
			@settings.new_block(blockname)
			@settings.block(blockname).copy_from!(@settings.block(:calculation), :keep)
			settings = @settings.block(blockname)
			# Modify settings
			settings[:demon_plane_ratio] = @scan_ratios[i]
			settings[:demon_no_scf] = true
			settings[:demon_fragment_guess] = true
			calculation = Calculation.new(@name + "_" + name, settings, @geometry)
			calculation.prepare([:energy, :atomic_charges])
			@calculations << calculation
		}
	end

	def run_protocol(queue)
		# Submit calculations
		@calculations.each{|calculation|
			calculation.send_to_queue(queue)
		}
		# Print the scan
		puts "Plane ratio scan:"
		@scan_populations = []
		@scan_ratios.each_index{|i|
			#!# Temporary solution: the charge is found on the first atom of the selection
			atom1 = @geometry.atomlist_from_selection(@calculations[i].settings[:demon_plane_axis_a])[0]
			@scan_populations[i] = @calculations[i].results.atomic_charges[atom1]
			printf("%8.3f %10.5f\n", @scan_ratios[i], @scan_populations[i])

		}
		# Check if there is sign change, save the interval
		ia = ib = nil
		(@scan_populations.size - 1).times{|i|
			if @scan_populations[i] * @scan_populations[i+1] < 0
				ia = i
				ib = i+1
			end
		}
		if ia.nil?
			Cuby::log.puts "\nNeutral fragments can not be found within the scanned interval\n\n"
			exit
		end

		# Interpolate the zero crossing
		xa = @scan_ratios[ia]
		xb = @scan_ratios[ib]
		ya = @scan_populations[ia].abs
		yb = @scan_populations[ib].abs
		@optimal_ratio = xa + (ya / (ya+yb)) * (xb - xa)
		Cuby::log.puts "Interpolated partitioning ratio: #{'%.4f' % @optimal_ratio}"
		Cuby::log.puts_v_only(:normal,"")

		# Prepare free calculation
		name = 'free'
		blockname = ("calculation_" + name).to_sym
		@settings.new_block(blockname)
		@settings.block(blockname).copy_from!(@settings.block(:calculation), :keep)
		settings = @settings.block(blockname)
		# Modify settings
		settings[:demon_plane_ratio] = @optimal_ratio
		@calculation_free = Calculation.new(@name + "_" + name, settings, @geometry)
		@calculation_free.prepare([:energy, :atomic_charges])
		@calculation_free.send_to_queue(queue)

		# Prepare constrained DFT calculation
		name = 'cdft'
		blockname = ("calculation_" + name).to_sym
		@settings.new_block(blockname)
		@settings.block(blockname).copy_from!(@settings.block(:calculation), :keep)
		settings = @settings.block(blockname)
		# Modify settings
		settings[:demon_plane_ratio] = @optimal_ratio
		settings[:demon_constrained_dft] = :spatial
		@calculation_cdft = Calculation.new(@name + "_" + name, settings, @geometry)
		@calculation_cdft.prepare([:energy, :atomic_charges])
		@calculation_cdft.send_to_queue(queue)
	end

	def print_protocol
		# Print results
		results.print_energy
		Cuby::log.puts_v_only(:normal,"")
		results.print_energy_decomposition if @settings[:print].include?(:energy_decomposition)
		results.print_dipole if @settings[:print].include?(:dipole)
	end

	def cleanup_protocol
		@calculation_free.cleanup
		#@calculation_cdft.cleanup
		@calculations.each{|calculation|
			calculation.cleanup
		}
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			# Get results
			@results = Results.new
			@results.energy = @calculation_free.results.energy - @calculation_cdft.results.energy
		end
		return @results
	end
end
