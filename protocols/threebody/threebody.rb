
module ProtocolThreebody
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :composite_calculation
	PROTOCOL_DESCRIPTION = "Calculation of threebody energy in non-covalent trimers"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:molecule_a, :optional, "Definition of the first monomer (selection, charge, multiplicity)"],
		InputBlock[:molecule_b, :optional, "Definition of the second monomer (selection, charge, multiplicity)"],
		InputBlock[:molecule_c, :optional, "Definition of the third monomer (selection, charge, multiplicity)"]
	]
	#=======================================================================

	# No results yet
	@results = nil

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Threebody effects calculation')
	end

	def prepare_protocol
		# Read geometry
		load_geometry

		# Fix settings
		if @settings.block(:molecule_a) && @settings.block(:molecule_b) && @settings.block(:molecule_c)
			# Both set
		elsif !@settings.block(:molecule_a) && !@settings.block(:molecule_b) && !@settings.block(:molecule_c)
			# None set - automated setup
			@settings[:molecule_a, :selection] = "auto"
			@settings[:molecule_b, :selection] = "auto"
			@settings[:molecule_c, :selection] = "auto"
			# Charge
			if !@settings.set?(:charge) || @settings[:charge] == 0
				Cuby::warning("Charges not defined for monomers, using default value \"0\"")
				@settings[:molecule_a, :charge] = 0
				@settings[:molecule_b, :charge] = 0
				@settings[:molecule_c, :charge] = 0
			else
				Cuby::error("Charge of monomers can not be automatically assigned")
			end
			# Multiplicity
			if @settings[:multiplicity] != 1
				Cuby::error("Multiplicity of monomers can not be automatically assigned")
			end
		else
			Cuby::error("Either three monomer calculations must be defined, or none for automatic input")
		end

		# Auto selections
		@settings[:molecule_a, :selection] = "%molecule(1;3)" if @settings[:molecule_a, :selection] == "auto"
		@settings[:molecule_b, :selection] = "%molecule(2;3)" if @settings[:molecule_b, :selection] == "auto"
		@settings[:molecule_c, :selection] = "%molecule(3;3)" if @settings[:molecule_c, :selection] == "auto"

		# Calculation scheme
		trimer    = ["abc"]
		dimers    = ["abx", "axc", "xbc"]
		monomers  = ["axx", "xbx", "xxc"]

		case @settings[:bsse_correction_3b]
		when :none
			@scheme = trimer + dimers + monomers
		when :full
			@scheme = trimer + dimers + monomers
			@scheme.map!{|s| s.gsub("x","g")}
		end
		@weights = [
			1.0, # trimer
			-1.0, -1.0, -1.0, # dimers
			1.0, 1.0, 1.0 # monomers
		]

		# Monomer geometries and charges
		geometries = []
		geometries_ghost = []
		charges = []
		[:molecule_a, :molecule_b, :molecule_c].each{|mol_x|
			geo = @geometry.geometry_from_selection(@settings[mol_x, :selection])
			geo_g = geo.deep_copy
			geo_g.each{|atom| atom.properties[:ghost] = true}

			geometries << geo
			geometries_ghost << geo_g
			charges << @settings[mol_x, :charge]
		}

		# Block names
		@calculation_names = []
		@scheme.each{|calcname|
			@calculation_names << "calculation_#{calcname}".to_sym
		}

		# Build calculations
		@calculations = []
		@scheme.each{|calcname|
			# Create block
			settings = @settings.new_block("calculation_#{calcname}".to_sym)
			settings.copy_from!(@settings, :keep, [:molecule_a, :molecule_b, :molecule_c] + @calculation_names)
			geometry = Geometry.new
			settings[:charge] = 0
			calcname.size.times{|i|
				letter = calcname[i..i]
				if letter =~ /[a-c]/
					# Use the system
					geometry = geometry + geometries[i]
					settings[:charge] += charges[i]
				elsif letter == "g"
					# Use only basis
					geometry = geometry + geometries_ghost[i]
				elsif letter == "x"
					# Not present
				end
			}
			calculation = Calculation.new(@name + calcname, settings, geometry)
			@calculations << calculation
		}

		# Prepare calculations
		@calculations.each{|calculation|
			calculation.prepare([:energy])
		}
	end

	def run_protocol(queue)
		# Submit calculation
		@calculations.each{|calculation|
			calculation.send_to_queue(queue)
		}
	end

	def print_protocol
		# Subsystem results
		if @settings[:print].include?(:subsystem_results)
			Cuby::log.puts_v_only(:normal,"Susbsystem energies:")
			@scheme.each_index{|i|
				name = "Subsystem #{@scheme[i].gsub(/[g,x]/,"").upcase}"
				@calculations[i].results.print_energy(name)
				@calculations[i].results.print_energy_decomposition(name) if @settings[:print].include?(:energy_decomposition)

			}
			Cuby::log.puts_v_only(:normal,"")
		end

		# Print results
		results.print_energy("3-body energy")
		Cuby::log.puts_v_only(:normal,"")
		results.print_energy_decomposition("3-body energy") if @settings[:print].include?(:energy_decomposition)
	end

	def cleanup_protocol
		@calculations.each{|calculation|
			calculation.cleanup
		}
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			# Get results
			@results = Results.new

			@results.energy = 0.0
			# Add energies with their weights
			@calculations.each_index{|i|
				@results.energy += @weights[i] * @calculations[i].results.energy
			}
			# Energy components
			@calculations[0].results.energy_components.each_pair{|key, value|
				c = 0.0
				@calculations.each_index{|i|
					c += @weights[i] * @calculations[i].results.energy_components[key]
				}
				@results.energy_components[key] = c
			}
		end
		return @results
	end

	#=======================================================================
	# Optional protocol functionality
	#=======================================================================
	
	def single_line_result
		return sprintf("3-body energy: %20.6f kcal/mol", results.energy)
	end

	def single_number_result
		return results.energy
	end
end
