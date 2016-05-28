################################################################################
#
# Mixer interface
#
# Author: Jan Rezac
# Date created: 2013-07-17
# License: Cuby4 license
# Description: Composite interface for combining two calculations
# Status: Works
#
################################################################################

#===============================================================================
# Mixes two calculations with user-defined weights
#
# For weights wA and wB, the resulting energy is
# E(mixed) = wA * E(A) + wB * E(B),
# the same is applied to gradient.
#
#===============================================================================

module InterfaceMixer
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
		InputBlock[:calculation_a, :required, "Setup for the first calculation (A)"],
		InputBlock[:calculation_b, :required, "Setup for the second calculation (B)"],
		InputBlock[:calculation_common, :optional, "Setup shared by both calculations"]
	]
	#=======================================================================

	def prepare_interface
		# Check child settings
		Cuby::error("Input block calculation_a must be defined") unless @settings.has_block?(:calculation_a)
		Cuby::error("Input block calculation_b must be defined") unless @settings.has_block?(:calculation_b)

		# Child settings
		calc_settings_a = @settings.block(:calculation_a)
		calc_settings_b = @settings.block(:calculation_b)

		calc_settings_a.copy_from!(@settings.block(:calculation_common), :keep) if @settings.has_block?(:calculation_common)
		calc_settings_b.copy_from!(@settings.block(:calculation_common), :keep) if @settings.has_block?(:calculation_common)
		# Copy charge if it is set at the root level but not in childs
		calc_settings_a[:charge] = @settings[:charge] if @settings.set?(:charge) && ! calc_settings_a.set?(:charge)
		calc_settings_b[:charge] = @settings[:charge] if @settings.set?(:charge) && ! calc_settings_b.set?(:charge)
		calc_settings_a[:multiplicity] = @settings[:multiplicity] if @settings.set?(:multiplicity) && ! calc_settings_a.set?(:multiplicity)
		calc_settings_b[:multiplicity] = @settings[:multiplicity] if @settings.set?(:multiplicity) && ! calc_settings_b.set?(:multiplicity)

		# Build calculations
		@calc_a = Calculation.new(@name+'_A', calc_settings_a, @geometry)
		@calc_b = Calculation.new(@name+'_B', calc_settings_b, @geometry)

		# Prepare calculations
		@calc_a.prepare(@what)
		@calc_b.prepare(@what)

	end

	def queue_interface(queue)
		@calc_a.send_to_queue(queue)
		@calc_b.send_to_queue(queue)
	end

	def compose_interface
		results = Results.new

		results.energy = @calc_a.results.energy * @settings[:mixer_weight_a] + @calc_b.results.energy *  @settings[:mixer_weight_b]

		results.energy_components[:energy_a] = @calc_a.results.energy
		results.energy_components[:energy_b] = @calc_b.results.energy
		
		if @what.include?(:gradient)
			results.gradient = 
				@calc_a.results.gradient * @settings[:mixer_weight_a] + 
				@calc_b.results.gradient * @settings[:mixer_weight_b] 
		end

		if @what.include?(:hessian)
			results.hessian = 
				@calc_a.results.hessian * @settings[:mixer_weight_a] + 
				@calc_b.results.hessian * @settings[:mixer_weight_b] 
		end

		return results
	end

	def cleanup_interface
		@calc_a.cleanup
		@calc_b.cleanup
	end
end
