
module ProtocolInteraction
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :composite_calculation
	PROTOCOL_DESCRIPTION = "Interaction energy in non-covalent complexes"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:molecule_a, :optional, "Definition of the first monomer (selection, charge, multiplicity)"],
		InputBlock[:molecule_b, :optional, "Definition of the second monomer (selection, charge, multiplicity)"]
	]
	#=======================================================================

	# No results yet
	@results = nil

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Interaction energy calculation')
	end

	def prepare_protocol
		# Read geometry
		load_geometry

		# Fix settings
		if @settings.block(:molecule_a) && @settings.block(:molecule_b)
			# Both set
		elsif !@settings.block(:molecule_a) && !@settings.block(:molecule_b)
			# None set - automated setup
			@settings[:molecule_a, :selection] = "auto"
			@settings[:molecule_b, :selection] = "auto"
			# Charge
			if !@settings.set?(:charge) || @settings[:charge] == 0
				Cuby::warning("Charges not defined for monomers, using default value \"0\"")
				@settings[:molecule_a, :charge] = 0
				@settings[:molecule_b, :charge] = 0
			else
				Cuby::error("Charge of monomers can not be automatically assigned")
			end
			# Multiplicity
			if @settings[:multiplicity] != 1
				Cuby::error("Multiplicity of monomers can not be automatically assigned")
			end
		else
			Cuby::error("Two monomers must be defined in the input.")
		end

		# Auto selections
		@settings[:molecule_a, :selection] = "%molecule(1;2)" if @settings[:molecule_a, :selection] == "auto"
		@settings[:molecule_b, :selection] = "%molecule(2;2)" if @settings[:molecule_b, :selection] == "auto"

		# Create settings for dimer
		@settings.new_block(:molecule_ab) unless @settings.has_block?(:molecule_ab)

		# Copy everything from root level to system settings, keeping what is set, skipping molecule_* blocks
		@settings.block(:molecule_a).copy_from!(@settings, :keep, [:molecule_a, :molecule_b, :molecule_ab])
		@settings.block(:molecule_b).copy_from!(@settings, :keep, [:molecule_a, :molecule_b, :molecule_ab])
		@settings.block(:molecule_ab).copy_from!(@settings, :keep, [:molecule_a, :molecule_b, :molecule_ab])

		# Total charge
		new_charge = @settings[:molecule_a, :charge] + @settings[:molecule_b, :charge]
		if @settings.set?(:molecule_ab, :charge)
			Cuby::warning("Total charge was different from sum of monomer charges") unless @settings[:molecule_ab, :charge] == new_charge
		end
		@settings[:molecule_ab, :charge] = new_charge

		# Monomer geometries
		list_a = @geometry.atomlist_from_selection(@settings[:molecule_a, :selection])
		list_b = @geometry.atomlist_from_selection(@settings[:molecule_b, :selection])
		if GeometrySelections.overlap?(list_a, list_b)
			Cuby::error("Monomer selections do overlap")
		end
		geometry_a = @geometry.geometry_from_list(list_a)
		geometry_b = @geometry.geometry_from_list(list_b)
		geometry_ab = geometry_a + geometry_b
		if geometry_ab.size != @geometry.size && @settings[:interaction_check_size]
			Cuby::warning("Some atoms from the input geometry are not used in the interaction calculation")
		end

		# BSSE correction: add ghost atoms
		if @settings[:bsse_correction]
			geometry_ag = geometry_a.deep_copy
			geometry_bg = geometry_b.deep_copy

			geometry_ag.each{|atom| atom.properties[:ghost] = true}
			geometry_bg.each{|atom| atom.properties[:ghost] = true}

			geometry_a = geometry_a + geometry_bg
			geometry_b = geometry_b + geometry_ag
		end

		# Create calculation
		@calculation_ab = Calculation.new(@name + "_ab", @settings.block(:molecule_ab), geometry_ab)
		@calculation_a =  Calculation.new(@name + "_a",  @settings.block(:molecule_a),  geometry_a)
		@calculation_b =  Calculation.new(@name + "_b",  @settings.block(:molecule_b),  geometry_b)
		# Prepare calculation
		@calculation_ab.prepare([:energy])
		@calculation_a.prepare([:energy])
		@calculation_b.prepare([:energy])
	end

	def run_protocol(queue)
		# Submit calculation
		@calculation_ab.send_to_queue(queue)
		@calculation_a.send_to_queue(queue)
		@calculation_b.send_to_queue(queue)
	end

	def print_protocol
		# Subsystem energies
		if @settings[:print].include?(:subsystem_results)
			Cuby::log.puts_v_only(:normal,"Susbsystem energies:")
			{"Dimer AB" => @calculation_ab, "Monomer A" => @calculation_a, "Monomer B" => @calculation_b}.each_pair {|name, calculation|
				calculation.results.print_energy(name)
				calculation.results.print_energy_decomposition(name) if @settings[:print].include?(:energy_decomposition)

			}
			Cuby::log.puts_v_only(:normal,"")
		end

		# Print results
		results.print_energy("Interaction energy")
		Cuby::log.puts_v_only(:normal,"")
		results.print_energy_decomposition("Interaction energy") if @settings[:print].include?(:energy_decomposition)
	end

	def cleanup_protocol
		@calculation_ab.cleanup
		@calculation_a.cleanup
		@calculation_b.cleanup
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			# Get results
			@results = Results.new
			@results.energy = @calculation_ab.results.energy - @calculation_a.results.energy - @calculation_b.results.energy
			# Energy components
			@calculation_ab.results.energy_components.each_pair{|key, value|
				@results.energy_components[key] = value - @calculation_a.results.energy_components[key] -  @calculation_b.results.energy_components[key]
			}
		end
		return @results
	end

	#=======================================================================
	# Optional protocol functionality
	#=======================================================================
	
	def single_line_result
		return sprintf("Interaction energy: %20.6f kcal/mol", results.energy)
	end

	def single_number_result
		return results.energy
	end
end
