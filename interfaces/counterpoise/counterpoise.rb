################################################################################
#
# Interface Counterpoise
#
# Author: Jan Rezac
# Date created: 2013-02-05
# License: Cuby4 license
# Description: Composite interface implementing the counterpoise correction for BSSE
# Status: Works
#
################################################################################

#===============================================================================
# Counterpoise correction of energy and its derivatives
#===============================================================================

module InterfaceCounterpoise
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :composite
	CAPABILITIES = [:energy, :gradient, :hessian]
	MODIFIER = false
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation, :required, "Setup for the calculation"],
		InputBlock[:molecule_a, :optional, "Block defining monomer A (should contain keywords selection, charge, multiplicity)"],
		InputBlock[:molecule_b, :optional, "Block defining monomer B (should contain keywords selection, charge, multiplicity)"]
	]
	#=======================================================================

	def prepare_interface
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
		@settings.new_block(:molecule_ab)

		# Create settings for monomers with ghost atoms
		@settings.new_block(:molecule_ag)
		@settings.new_block(:molecule_bg)

		# Copy everything from root level to system settings, keeping what is set, skipping molecule_* blocks
		@settings.block(:molecule_ab).copy_from!(@settings.block(:calculation), :keep, [/^molecule_/])
		@settings.block(:molecule_a).copy_from!(@settings.block(:calculation), :keep, [/^molecule_/])
		@settings.block(:molecule_b).copy_from!(@settings.block(:calculation), :keep, [/^molecule_/])
		@settings.block(:molecule_ag).copy_from!(@settings.block(:molecule_a), :keep)
		@settings.block(:molecule_bg).copy_from!(@settings.block(:molecule_b), :keep)

		# Total charge
		new_charge = @settings[:molecule_a, :charge] + @settings[:molecule_b, :charge]
		if @settings.set?(:molecule_ab, :charge)
			Cuby::warning("Total charge was different from sum of monomer charges") unless @settings[:molecule_ab, :charge] == new_charge
		end
		@settings[:molecule_ab, :charge] = new_charge

		# Monomer geometries
		@list_a = @geometry.atomlist_from_selection(@settings[:molecule_a, :selection])
		@list_b = @geometry.atomlist_from_selection(@settings[:molecule_b, :selection])
		if GeometrySelections.overlap?(@list_a, @list_b)
			Cuby::error("Monomer selections do overlap")
		end
		@geometry_a = @geometry.geometry_from_list(@list_a).deep_copy
		@geometry_b = @geometry.geometry_from_list(@list_b).deep_copy
		unless @geometry_a.size + @geometry_b.size == @geometry.size
			Cuby::error("Some atoms from the input geometry are not used in the counterpoise calculation")
		end

		# Geometries with ghost atoms
		@geometry_ag = ghosts_geo(@geometry_a)
		@geometry_bg = ghosts_geo(@geometry_b)

		# Create calculations
		@calculation_ab = Calculation.new(@name + "_ab", @settings.block(:molecule_ab), @geometry_a + @geometry_b)
		@calculation_a =  Calculation.new(@name + "_a",  @settings.block(:molecule_a),  @geometry_a)
		@calculation_b =  Calculation.new(@name + "_b",  @settings.block(:molecule_b),  @geometry_b)
		@calculation_ag =  Calculation.new(@name + "_ag",  @settings.block(:molecule_ag),  @geometry_a + @geometry_bg)
		@calculation_bg =  Calculation.new(@name + "_bg",  @settings.block(:molecule_bg),  @geometry_ag + @geometry_b)

		@calculations = [@calculation_ab, @calculation_ag, @calculation_bg, @calculation_a, @calculation_b]

		# Prepare calculations
		@calculations.each{|calc| calc.prepare(@what)}

	end

	def queue_interface(queue)
		# Update geometries
		update_child_geometries!

		# Send calculations to queue
		@calculations.each{|calc| calc.send_to_queue(queue)}
	end

	def compose_interface
		results = Results.new

		results.energy = counterpoise_energy
		counterpoise_energy_components(results.energy_components)
		if @what.include?(:gradient)
			results.gradient = counterpoise_gradient
		end
		if @what.include?(:hessian)
			results.hessian = counterpoise_hessian
		end
		return results
	end

	def cleanup_interface
		@calculations.each{|calc| calc.cleanup}
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def ghosts_geo(ghostgeo)
		ghosts = ghostgeo.deep_copy
		ghosts.each{|atom| atom.properties[:ghost] = true}
		return ghosts
	end

	def update_child_geometries!
		@list_a.each_index{|i|
			@geometry_a[i].set_coord(@geometry[@list_a[i]])
			@geometry_ag[i].set_coord(@geometry[@list_a[i]])
		}
		@list_b.each_index{|i|
			@geometry_b[i].set_coord(@geometry[@list_b[i]])
			@geometry_bg[i].set_coord(@geometry[@list_b[i]])
		}
	end

	def counterpoise_energy
		return @calculation_ab.results.energy -
			@calculation_ag.results.energy -
			@calculation_bg.results.energy +
			@calculation_a.results.energy +
			@calculation_b.results.energy
	end

	def counterpoise_energy_components(components)
		# Calculate counterpoise-corrected energy components
		@calculation_ab.results.energy_components.each_key{|key|
			components[key] = @calculation_ab.results.energy_components[key] -
				@calculation_ag.results.energy_components[key] -
				@calculation_bg.results.energy_components[key] +
				@calculation_a.results.energy_components[key] +
				@calculation_b.results.energy_components[key]
		}

		# Add new components
		components[:energy_nocp] = @calculation_ab.results.energy
		components[:interaction_energy_nocp] = @calculation_ab.results.energy - @calculation_a.results.energy - @calculation_b.results.energy
		components[:interaction_energy_cp] = @calculation_ab.results.energy - @calculation_ag.results.energy - @calculation_bg.results.energy

		return nil
	end

	def counterpoise_gradient
		grad_a = @calculation_a.results.gradient.append(Gradient.zero(@geometry_b.size))
		grad_b = Gradient.zero(@geometry_a.size).append(@calculation_b.results.gradient)

		composite = @calculation_ab.results.gradient -
			@calculation_ag.results.gradient - 
			@calculation_bg.results.gradient +
			grad_a + grad_b

		# The gradient must be reordered to match the original geometry
		reordered =  Gradient.new
		list = @list_a + @list_b
		list.each_index{|i|
			reordered[list[i]] = composite[i]
		}

		return reordered
	end

	def counterpoise_hessian
		#!# Correct for different ordering
		h_subsystems = Hessian.zero(@calculation_ab.results.hessian.n)

		h_subsystems.paste!(0, 0, @calculation_a.results.hessian)
		h_subsystems.paste!(@calculation_a.results.hessian.m, @calculation_a.results.hessian.m, @calculation_b.results.hessian)

		composite = @calculation_ab.results.hessian -
			@calculation_ag.results.hessian -
			@calculation_bg.results.hessian +
			h_subsystems

		# The Hessian must be reordered to match the original geometry
		reordered = Hessian.zero(@calculation_ab.results.hessian.n)
		list = @list_a + @list_b
		list.each_index{|i|
			list.each_index{|j|
				element = composite.submatrix(i*3, j*3, 3, 3)
				reordered.paste!(list[i]*3, list[j]*3, element)
			}
		}

		return reordered
	end
end
