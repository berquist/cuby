################################################################################
#
# Extrapolation to CBS interface
#
# Author: Jan Rezac
# Date created: 2013-04-14
# License: Cuby4 license
# Description: Composite interface implementing the MP2 and CCSD(T)/CBS extrapolation
# Status: Works
#
################################################################################

#===============================================================================
# Implements MP2/CBS extrapolation + higher-order correction, eg. composite
# CCSD(T)/CBS
#===============================================================================

module InterfaceExtrapolateCbs
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :composite
	CAPABILITIES = [:energy, :gradient, :hessian, :ghost_atoms]
	MODIFIER = false
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation_mp2_small, :required, "Input for the MP2 calculation in the smaller basis set"],
		InputBlock[:calculation_mp2_large, :required, "Input for the MP2 calculation in the larger basis set"],
		InputBlock[:calculation_corr,      :optional, "Input for the higher-order correction, e.g. CCSD(T)"],
		InputBlock[:calculation_hf_small,  :optional, "For gradient calculations, a separate HF calculation in the smaller basis set is needed"],
		InputBlock[:calculation_hf_large,  :optional, "For gradient calculations, a separate HF calculation in the larger basis set is needed"],
		InputBlock[:calculation_corr_mp2,  :optional, "For gradient calculations, a separate MP2 calculation in the basis set of the higher-order correction"],
		InputBlock[:calculation_common,    :optional, "Common setup, applied to all the calculations"]
	]
	#=======================================================================

	# Child calculations and their abbreviations
	CHILD_NAMES = {
		:calculation_corr	=> :Corr,
		:calculation_corr_mp2	=> :CorrMP2,
		:calculation_mp2_small	=> :MP2s,
		:calculation_mp2_large	=> :MP2l,
		:calculation_hf_small	=> :HFs,
		:calculation_hf_large	=> :HFl
	}

	def prepare_interface
		if @what.include?(:gradient) || @what.include?(:hessian)
			case @settings[:extrapolate_cbs_grad_mode]
			when :use_mp2_gradient_components
				@calc_names = [:calculation_mp2_small, :calculation_mp2_large]
				if @what.include?(:hessian)
					Cuby::error("Hessian calculation using extrapolate_cbs interface possible only\nwhen 'extrapolate_cbs_grad_mode' is set to 'separate_calculations'")
				end
			when :separate_calculations
				@calc_names = [:calculation_mp2_small, :calculation_mp2_large, :calculation_hf_small, :calculation_hf_large]
			end
			if @settings[:extrapolate_cbs_correction]
				@calc_names << :calculation_corr
				@calc_names << :calculation_corr_mp2
			end
		else
			@calc_names = [:calculation_mp2_small, :calculation_mp2_large]
			if @settings[:extrapolate_cbs_correction]
				@calc_names << :calculation_corr
			end
		end

		# The required blocks
		@calc_names.each{|blockname|
			unless @settings.has_block?(blockname)
				Cuby::error("Block \"#{blockname}\" not found in the input")
			end
		}

		# Copy from calculation_common block
		@calc_names.each{|blockname|
			@settings.block(blockname).copy_from!(@settings.block(:calculation_common), :keep) if @settings.has_block?(:calculation_common)
			# Set charge to parent if not defined explicitly
			@settings.block(blockname)[:charge] = @settings[:charge] unless @settings.block(blockname).set?(:charge)
			@settings.block(blockname)[:multiplicity] = @settings[:multiplicity] unless @settings.block(blockname).set?(:multiplicity)
		}


		# Build calculations
		@calculations = {}
		@calc_names.each{|blockname|
			name = CHILD_NAMES[blockname]
			@calculations[name] = Calculation.new(@name+"_"+name.to_s, @settings.block(blockname), @geometry)

		}

		# Prepare calculations
		@calculations.each_value{|calc| calc.prepare(@what)}
	end

	def queue_interface(queue)
		@calculations.each_value{|calc| calc.send_to_queue(queue)}
	end

	def compose_interface
		results = Results.new

		ehf = @calculations[:MP2l].results.energy_components[:scf_energy]

		emp2s = @calculations[:MP2s].results.energy_components[:correlation_mp2]
		emp2l = @calculations[:MP2l].results.energy_components[:correlation_mp2]


		bs = @settings[:calculation_mp2_small, :basisset_zeta]
		bl = @settings[:calculation_mp2_large, :basisset_zeta]

		# SCF energy - taken from the larger MP2 calculation
		results.energy_components[:scf_energy] = ehf
		# MP2 correlation energy - Helgaker extrapolation
		results.energy_components[:correlation_mp2_cbs] = (emp2l * bl**3 - emp2s * bs**3)/(bl**3 - bs**3)

		delta_corr_weight = @settings[:extrapolate_cbs_corr_weight]

		if @settings[:extrapolate_cbs_correction]
			# Calculate CCSD(T) correction
			if @calc_names.include?(:calculation_corr_mp2)
				# dCC from two separate calculations
				delta_corr =  @calculations[:Corr].results.energy - @calculations[:CorrMP2].results.energy
			else
				# dCC from energy components
				eccmp2 = @calculations[:Corr].results.energy_components[:scf_energy] + @calculations[:Corr].results.energy_components[:correlation_mp2]
				ecc =    @calculations[:Corr].results.energy
				delta_corr =  ecc - eccmp2
			end
			# CCSD(T)/CBS energy
			results.energy_components[:"correlation_composite_cbs"] = results.energy_components[:correlation_mp2_cbs] + delta_corr * delta_corr_weight
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:"correlation_composite_cbs"]
		else
			# Only MP2/CBS
			results.energy = results.energy_components[:scf_energy] + results.energy_components[:"correlation_mp2_cbs"]
		end


		if @what.include?(:gradient)
			# MP2 correlation gradients and SCF gradient
			case @settings[:extrapolate_cbs_grad_mode]
			when :use_mp2_gradient_components
				# Check gradient components
				if @calculations[:MP2l].results.gradient_components.empty? || @calculations[:MP2s].results.gradient_components.empty?
					Cuby::error "Gradient components not available, use interface that supports this feature\nand switch it on in the input."
				end
				# Components available
				gmp2l = @calculations[:MP2l].results.gradient_components[:correlation_gradient]
				gmp2s = @calculations[:MP2s].results.gradient_components[:correlation_gradient]
				gscf = @calculations[:MP2l].results.gradient_components[:scf_gradient]
			when :separate_calculations
				# Correlation gradients from separate calculations
				gmp2l = @calculations[:MP2l].results.gradient - @calculations[:HFl].results.gradient
				gmp2s = @calculations[:MP2s].results.gradient - @calculations[:HFs].results.gradient
				gscf = @calculations[:HFl].results.gradient
			end
			# Compose
			results.gradient = gscf + (gmp2l * bl**3 - gmp2s * bs**3) * (1.0/(bl**3 - bs**3))
			# dCCSD(T) gradient
			if @settings[:extrapolate_cbs_correction]
				g_delta_corr = @calculations[:Corr].results.gradient - @calculations[:CorrMP2].results.gradient
				results.gradient = results.gradient + g_delta_corr * delta_corr_weight
			end
		end

		if @what.include?(:hessian)
			case @settings[:extrapolate_cbs_grad_mode]
			when :separate_calculations
				# SCF hessian
				h_scf = @calculations[:HFl].results.hessian
				# MP2 correlation hessians
				h_mp2l = @calculations[:MP2l].results.hessian - @calculations[:HFl].results.hessian
				h_mp2s = @calculations[:MP2s].results.hessian - @calculations[:HFs].results.hessian
			else
				raise "Weird error, this code shouldn't be accessible for now"
			end
			results.hessian = h_scf + ((h_mp2l * bl**3) - (h_mp2s * bs**3)) * (1.0/(bl**3 - bs**3))
			# dCCSD(T) hessian
			if @settings[:extrapolate_cbs_correction]
				h_delta_corr = @calculations[:Corr].results.hessian - @calculations[:CorrMP2].results.hessian
				results.hessian = results.hessian + h_delta_corr * delta_corr_weight
			end
		end

		return results
	end

	def cleanup_interface
		@calculations.each_value{|calc| calc.cleanup}
	end
end
