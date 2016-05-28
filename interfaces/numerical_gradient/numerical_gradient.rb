################################################################################
#
# Numerical gradient interface
#
# Author: Jan Rezac
# Date created: 2012-04-11
# License: Cuby4 license
# Description: Numerical gradient calculation
# Status: Works
#
################################################################################

#===============================================================================
# Numerical gradient calculation via finite differences
#===============================================================================

module InterfaceNumericalGradient
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :composite
	CAPABILITIES = [:energy, :gradient]
	MODIFIER = false
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation, :required, "Setup of the actual calculation"]
	]
	#=======================================================================

	def prepare_interface

		case @settings[:numerical_gradient_type]
		when :single_sided
			@displacements = [
				Coordinate[1,0,0],
				Coordinate[0,1,0],
				Coordinate[0,0,1]
			]
		when :central
			@displacements = [
				Coordinate[ 1, 0, 0],
				Coordinate[-1, 0, 0],
				Coordinate[ 0, 1, 0],
				Coordinate[ 0,-1, 0],
				Coordinate[ 0, 0, 1],
				Coordinate[ 0, 0,-1]
			]
		when :tetrahedral
			@displacements = [
				Coordinate[ 1, 1, 1],
				Coordinate[ 1,-1,-1],
				Coordinate[-1, 1,-1],
				Coordinate[-1,-1, 1]
			]
		end

		# Child settings
		calc_settings = @settings.block(:calculation)

		child_what = [:energy]
		child_what << :atomic_charges if @settings[:numerical_gradient_atomic_charges]

		# Build child calculations

		# Energy is calculated on the original geometry
		# Used also as a reference for some differentiation schemes
		@energy_calc = Calculation.new(@name+'_0_000', calc_settings, @geometry)
		@energy_calc.prepare(child_what)

		# Build displaced geometries
		if @what.include?(:gradient)
			@grad_calcs = []
			@geometry.each_index{|i|
				@grad_calcs[i] = []

				@displacements.each_index{|di|
					d_vec = @displacements[di]
					geo = @geometry.deep_copy
					geo[i].plus!(d_vec * @settings[:numerical_gradient_step])
					crdstring = ""
					3.times{|c| crdstring += {1.0 => '+', -1.0 => '-', 0.0 => '0'}[d_vec[c]] }
					calc = Calculation.new(@name + "_#{i}_" + crdstring, calc_settings, geo)
					calc.prepare(child_what)
					@grad_calcs[i] << calc
				}
			}
		end
	end

	def queue_interface(queue)
		# queue energy calculation
		@energy_calc.send_to_queue(queue)

		if @what.include?(:gradient)
			# update geometries
			@geometry.each_index{|i|
				@displacements.each_index{|di|
					d_vec = @displacements[di]
					@grad_calcs[i][di].geometry.copy_coordinates_from(@geometry)
					@grad_calcs[i][di].geometry[i].plus!(d_vec * @settings[:numerical_gradient_step])
					geo = @geometry.deep_copy
					geo[i].plus!(d_vec * @settings[:numerical_gradient_step])
				}
			}

			# queue displaced calculations
			@grad_calcs.each{|a|
				a.each{|calc|
					calc.send_to_queue(queue)
				}
			}
		end
	end

	def compose_interface
		results = Results.new
		if @what.include?(:energy)
			results.energy = @energy_calc.results.energy
		end


		if @what.include?(:gradient)
			results.gradient = Gradient.new

			case @settings[:numerical_gradient_type]
			when :single_sided
				divider = 1.0
			when :central
				divider = 2.0
			when :tetrahedral
				divider = 4.0
			end

			@grad_calcs.each_index{|i|
				g = Coordinate.new
				@displacements.each_index{|di|
					d_vec = @displacements[di]
					ep = @grad_calcs[i][di].results.energy
					g += d_vec * (ep - @energy_calc.results.energy) / @settings[:numerical_gradient_step] / divider
				}
				results.gradient << g
			}
		end

		# Gradient of atomic charges
		if @settings[:numerical_gradient_atomic_charges]
			charge_gradient = AtomicChargesGradient.new(@geometry.size)
			@grad_calcs.each_index{|i| # atom index i
				@geometry.each_index{|j| # charge index j
					g = Coordinate.new
					@displacements.each_index{|di|
						d_vec = @displacements[di]
						ep = @grad_calcs[i][di].results.atomic_charges[j]
						g += d_vec * (ep - @energy_calc.results.atomic_charges[j]) / @settings[:numerical_gradient_step] / divider
					}
					charge_gradient[j][i] = g
					
				} 

			}

			results.atomic_charges_gradient = charge_gradient
		end

		return results
	end

	def cleanup_interface
		@energy_calc.cleanup

		if @what.include?(:gradient)
			@grad_calcs.each{|a|
				a.each{|calc|
					calc.cleanup
				}
			}
		end
	end
end
