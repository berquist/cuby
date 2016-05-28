require "protocols/heat_of_formation/classes/atom_thermochemistry.rb"

module ProtocolHeatOfFormation
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :composite_calculation
	PROTOCOL_DESCRIPTION = "Calculation of heat of formation of a molecule"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation_vib, :required, "Setup for the calculation of vibrations"],
		InputBlock[:calculation_ene, :required, "Setup for the calculation of energy"]
	]
	#=======================================================================
	
	@@element_data = nil
	ELEMENT_DATA_FILE = "data/atom_data_JANAF.yaml"

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Heat of formation')
	end

	def prepare_protocol
		# Input check
		unless @settings.has_block?(:calculation_vib)
			Cuby::error("Input for a Heat of formation job must contain subsection named\n'calculation_vib' which contains the setup for vibrational analysis")
		end
		unless @settings.has_block?(:calculation_ene)
			Cuby::error("Input for a Heat of formation job must contain subsection named\n'calculation_ene' which contains the setup for energy calculation")
		end

		# Read geometry
		load_geometry

		# Get elements in the system
		@elements = @geometry.count_elements

		# Load the table of element data
		unless @@element_data
			f = File.open(protocol_dir + '/' + ELEMENT_DATA_FILE)
			@@element_data = YAML.load(f)
			f.close
		end

		# Check if all elements are there
		@elements.each_key{|ele|
			Cuby::error("Thermochemical data not available for element #{ele}") unless @@element_data[ele]
		}
		
		# Get element multiplicities
		element_multiplicity = {}
		@elements.each_key{|ele|
			configuration = PeriodicTable.electron_configuration(ele)
			unpaired = 0
			configuration.orbitals.each_index{|i|
				orb = configuration.orbitals[i]
				pop = configuration.occupations[i]
				case orb
				when "s"
					pop = 0 if pop == 2
				when "p"
					pop -= 3 if pop > 3
				when "d"
					pop -= 5 if pop > 5
				when "f"
					pop -= 7 if pop > 7
				end
				unpaired += pop

			}
			element_multiplicity[ele] = unpaired + 1
		}

		# Build vibrational calculation
		settings = @settings.block(:calculation_vib)
		settings[:job] = :frequencies
		settings[:geometry] = @settings[:geometry]
		@job_vib = Job.new(@name + "_vib", settings)
		@job_vib.prepare

		# Build energy calculations
		@calculations_element = {}
		@calculations = []

		# First, the elements
		@elements.each_key{|ele|
			blockname = ("calculation_" + ele.to_s).to_sym
			@settings.new_block(blockname)
			@settings.block(blockname).copy_from!(@settings.block(:calculation_ene), :keep)
			@settings.block(blockname).copy_from!(@settings.block(:calculation_ene_element), :keep) if @settings.has_block?(:calculation_ene_element)
			settings = @settings.block(blockname)
			settings[:charge] = 0
			settings[:multiplicity] = element_multiplicity[ele]
			geometry = Geometry.new
			geometry << Atom.new(ele)
			calculation = Calculation.new(@name + "_" + ele.to_s, settings, geometry)
			calculation.prepare([:energy])
			@calculations << calculation
			@calculations_element[ele] = calculation
		}

		# Second, the whole system
		if @settings[:hof_energy_from_vibrations]
			@calculation_ene = @job_vib
		else
			settings = @settings.block(:calculation_ene)
			@calculation_ene = Calculation.new(@name + "_mol", settings, @geometry)
			@calculation_ene.prepare([:energy])
			@calculations << @calculation_ene
		end
	end

	def run_protocol(queue)
		# Submit vibrational calculation
		@job_vib.run

		# Submit energy calculations
		@calculations.each{|calculation|
			calculation.send_to_queue(queue)
		}
	end

	def print_protocol
		# Print results
		Cuby::log.puts_v_only(:normal,"--------------------------------------------------------------------------------")
		Cuby::log.puts_v_only(:normal,"Vibrational analysis")
		Cuby::log.puts_v_only(:normal,"-------------------------------------------------------------------------------")
		@job_vib.print_protocol
		Cuby::log.puts_v_only(:normal,"")

		Cuby::log.puts_v_only(:normal,"--------------------------------------------------------------------------------")
		Cuby::log.puts_v_only(:normal,"Molecule energy")
		Cuby::log.puts_v_only(:normal,"--------------------------------------------------------------------------------")
		@calculation_ene.results.print_energy
		Cuby::log.puts_v_only(:normal,"")

		Cuby::log.puts_v_only(:normal,"--------------------------------------------------------------------------------")
		Cuby::log.puts_v_only(:normal,"Element energies")
		Cuby::log.puts_v_only(:normal,"--------------------------------------------------------------------------------")
		@calculations_element.each_pair{|element, calculation|
			Cuby::log.puts_v_only(:normal, sprintf("%4s %16.6f kcal/mol", element.to_s, calculation.results.energy))
		}
		Cuby::log.puts_v_only(:normal,"")

		Cuby::log.puts_v_only(:normal,"--------------------------------------------------------------------------------")
		Cuby::log.puts_v_only(:normal,"Thermochemistry")
		Cuby::log.puts_v_only(:normal,"--------------------------------------------------------------------------------")
		results.print_energy("Heat of formation")
		
	end

	def cleanup_protocol
		@job_vib.cleanup
		@calculations.each{|calculation|
			calculation.cleanup
		}
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			# Get results
			@results = Results.new
			@results.energy = heat_of_formation
		end
		return @results
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def heat_of_formation
		# Atomization energy
		e_at = -1.0 * (@calculation_ene.results.energy + @job_vib.vibrations.zpve) # - ZPVE-corrected energy of the molecule
		@calculations_element.each_pair{|element, calculation| # + atomic energies
			e_at += @elements[element] * calculation.results.energy
		}
		@results.energy_components[:atomization_energy] = e_at

		# Heat of formation at 0K, deltaH_f(0K)
		hf_0 = 0.0
		@elements.each_pair{|element, count| # Sum of atomic corrections at 0K
			hf_0 += count * @@element_data[element].h0
		}
		hf_0 -= e_at
		@results.energy_components[:heat_of_formation_0K] = hf_0

		# Heat of formation at 298.15K
		hf_298 = hf_0
		hf_298 += (
			@job_vib.thermodynamics.energy_translation + 
			@job_vib.thermodynamics.energy_rotation + 
			@job_vib.thermodynamics.energy_vibration + 
			GAS_CONSTANT_SI * @job_vib.thermodynamics.temperature
		) * J2KCAL - @job_vib.vibrations.zpve # Thermal correction to H
		@elements.each_pair{|element, count| # Thermal correction for elements
			hf_298 -= count * @@element_data[element].h298_corr
		}

		return hf_298
	end
end
