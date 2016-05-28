################################################################################
#
# Numerical hessian interface
#
# Author: Jan Rezac
# Date created: 2013-03-18
# License: Cuby4 license
# Description: Numerical hessian calculation
# Status: Works
#
################################################################################

require "classes/math/polynomial_least_squares.rb"

#===============================================================================
# Numerical hessian calculation via finite differences
#===============================================================================

module InterfaceNumericalHessian
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
		InputBlock[:calculation, :required, "Setup of the actual calculation"]
	]
	#=======================================================================

	def prepare_interface
		return prepare_e_only if @settings[:numerical_hessian_type] == :energy_only

		case @settings[:numerical_hessian_type]
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
		else
			raise "Unknown numerical hessian type"
		end

		# Child settings
		calc_settings = @settings.block(:calculation)

		# Build child calculations

		# Energy is calculated on the original geometry
		# Used also as a reference for some differentiation schemes
		@energy_calc = Calculation.new(@name+'_0_000', calc_settings, @geometry)
		@energy_calc.prepare([:energy, :gradient])

		# Build displaced geometries
		if @what.include?(:hessian)
			@grad_calcs = []
			@geometry.each_index{|i|
				@grad_calcs[i] = []

				@displacements.each_index{|di|
					d_vec = @displacements[di]
					geo = @geometry.deep_copy
					geo[i].plus!(d_vec * @settings[:numerical_hessian_step])
					crdstring = ""
					3.times{|c| crdstring += {1.0 => '+', -1.0 => '-', 0.0 => '0'}[d_vec[c]] }
					calc = Calculation.new(@name + "_#{i}_" + crdstring, calc_settings, geo)
					calc.prepare([:energy, :gradient])
					@grad_calcs[i] << calc
				}
			}
		end
	end

	def queue_interface(queue)
		return queue_e_only(queue) if @settings[:numerical_hessian_type] == :energy_only

		# queue energy calculation
		@energy_calc.send_to_queue(queue)

		if @what.include?(:hessian)
			# update geometries
			@geometry.each_index{|i|
				@displacements.each_index{|di|
					d_vec = @displacements[di]
					@grad_calcs[i][di].geometry.copy_coordinates_from(@geometry)
					@grad_calcs[i][di].geometry[i].plus!(d_vec * @settings[:numerical_hessian_step])
					geo = @geometry.deep_copy
					geo[i].plus!(d_vec * @settings[:numerical_hessian_step])
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
		return compose_e_only if @settings[:numerical_hessian_type] == :energy_only

		results = Results.new
		if @what.include?(:energy)
			results.energy = @energy_calc.results.energy
		end

		if @what.include?(:gradient)
			results.gradient = @energy_calc.results.gradient
		end

		if @what.include?(:hessian)
			results.hessian = Hessian.zero(@geometry.size * 3)

			case @settings[:numerical_hessian_type]
			when :single_sided
				divider = 1.0
			when :central
				divider = 2.0
			when :tetrahedral
				divider = 4.0
			end

			g0 = @energy_calc.results.gradient

			@grad_calcs.each_index{|i|
				@displacements.each_index{|di|
					disp_v = Vector.zero(@geometry.size * 3)
					3.times{|k|
						disp_v[i*3 + k] = @displacements[di][k]
					}
					dgrad_v = (@grad_calcs[i][di].results.gradient - g0).to_vector

					results.hessian = results.hessian + ((disp_v.to_matrix * dgrad_v.to_matrix.transpose) * (1.0 / @settings[:numerical_hessian_step] / divider))
				}
			}

			# Symmetrize hessian
			results.hessian = (results.hessian + results.hessian.transpose) / 2

			# Dipole derivatives
			if @energy_calc.results.multipoles[:dipole]
				results.multipoles[:dipole] = @energy_calc.results.multipoles[:dipole]
				dipole_deriv = []
				results.multipoles[:dipole].derivative = dipole_deriv
				@grad_calcs.each_index{|i|
					#!# This has to be checked!
					3.times{|k|
						g = Coordinate.new
						@displacements.each_index{|di|
							d_vec = @displacements[di]
							dd = @grad_calcs[i][di].results.multipoles[:dipole][k] -
							     @energy_calc.results.multipoles[:dipole][k]
							g += d_vec * dd / @settings[:numerical_hessian_step] / divider
						}
						dipole_deriv << g
					}
				}
			end
		end

		return results
	end

	def cleanup_interface
		return cleanup_e_only if @settings[:numerical_hessian_type] == :energy_only

		@energy_calc.cleanup

		if @what.include?(:hessian)
			@grad_calcs.each{|a|
				a.each{|calc|
					calc.cleanup
				}
			}
		end
	end

	#========================================================================
	# Hessian from energy
	#========================================================================
	
	# 2nd derivative formulas from
	# http://www.holoborodko.com/pavel/numerical-methods/numerical-derivative/central-differences/
	
	def prepare_e_only
		# Child settings
		calc_settings = @settings.block(:calculation)

		# Build child calculations

		# Energy is calculated on the original geometry
		# This is also the zero for the differentiation
		@energy_calc = Calculation.new(@name+'_0_000', calc_settings, @geometry)
		@energy_calc.prepare([:energy])

		case @settings[:numerical_hessian_e_order]
		when 2
			@displacements = [[-1.0, -1.0],[1.0, 1.0],[-1.0, 1.0],[1.0, -1.0]]
			@displacements_g = [-1.0, 1.0]
		when 4
			@displacements = [
				[ 1.0, -2.0], [ 2.0, -1.0], [-2.0,  1.0], [-1.0, 2.0],
				[-1.0, -2.0], [-2.0, -1.0], [ 1.0,  2.0], [ 2.0, 1.0],
				[ 2.0, -2.0], [-2.0,  2.0], [-2.0, -2.0], [ 2.0, 2.0],
				[-1.0, -1.0], [ 1.0,  1.0], [ 1.0, -1.0], [-1.0, 1.0]
			]
			@displacements_g = [-2.0, -1.0, 1.0, 2.0]
		else
			Cuby::error("Value of the keyword numerical_hessian_e_order must be ither 2 or 4")
		end

		if @what.include?(:gradient)
			@displaced_calcs_g = []
			(@geometry.size * 3).times{|i|
				@displaced_calcs_g[i] = []
				@displacements_g.each_index{|d|
					geo = @geometry.deep_copy
					geo[i/3][i%3] += @displacements_g[d] * @settings[:numerical_hessian_step]
					calc = Calculation.new(@name + "_#{i}_#{d}", calc_settings, geo)
					calc.prepare([:energy])
					@displaced_calcs_g[i][d] = calc
				}
			}
		end

		if @what.include?(:hessian)
			# Displaced points indexing: coordinate a, coordinate b, displacements
			@displaced_calcs = []
			(@geometry.size * 3).times{|i|
				@displaced_calcs[i] = []
				(i+1).times{|j|
					@displaced_calcs[i][j] = []
					@displacements.each_index{|d|
						geo = @geometry.deep_copy
						geo[i/3][i%3] += @displacements[d][0] * @settings[:numerical_hessian_step]
						geo[j/3][j%3] += @displacements[d][1] * @settings[:numerical_hessian_step]
						calc = Calculation.new(@name + "_#{i}_#{j}_#{d}", calc_settings, geo)
						calc.prepare([:energy])
						@displaced_calcs[i][j][d] = calc
					}
				}
			}
		end
	end

	def queue_e_only(queue)
		# Queue the energy calculation
		@energy_calc.send_to_queue(queue)

		if @what.include?(:gradient)
			# Update geometries and queue displaced calculations
			(@geometry.size * 3).times{|i|
				@displacements_g.each_index{|d|
					geo = @displaced_calcs_g[i][d].geometry
					geo.copy_coordinates_from(@geometry)
					geo[i/3][i%3] += @displacements_g[d] * @settings[:numerical_hessian_step]
					@displaced_calcs_g[i][d].send_to_queue(queue)
				}
			}
		end

		if @what.include?(:hessian)
			# Update geometries and queue displaced calculations
			(@geometry.size * 3).times{|i|
				(i+1).times{|j|
					@displacements.each_index{|d|
						geo = @displaced_calcs[i][j][d].geometry
						geo.copy_coordinates_from(@geometry)
						geo[i/3][i%3] += @displacements[d][0] * @settings[:numerical_hessian_step]
						geo[j/3][j%3] += @displacements[d][1] * @settings[:numerical_hessian_step]
						@displaced_calcs[i][j][d].send_to_queue(queue)
					}
				}
			}
		end
	end

	def compose_e_only
		results = Results.new

		results.energy = @energy_calc.results.energy

		if @what.include?(:gradient)
			results.gradient = Gradient.zero(@geometry.size)
			(@geometry.size * 3).times{|i|
				y = []
				@displacements_g.each_index{|d|
					y[d] = @displaced_calcs_g[i][d].results.energy
				}
				case @settings[:numerical_hessian_e_order]
				when 2
					dfdx = (y[1] - y[0]) / 
					       (2.0 * @settings[:numerical_hessian_step])
				when 4
					dfdx = (y[0] - 8.0 * y[1] + 8.0 * y[2] - y[3]) /
					       (12.0 * @settings[:numerical_hessian_step])
				end
				results.gradient[i/3][i%3] = dfdx
			}
		end

		if @what.include?(:hessian)
			results.hessian = Hessian.zero(@geometry.size * 3)
			(@geometry.size * 3).times{|i|
				(i+1).times{|j|
					y = []
					@displacements.each_index{|d|
						y[d] = @displaced_calcs[i][j][d].results.energy
					}

					case @settings[:numerical_hessian_e_order]
					when 2
						d2dxy = 1.0/(4.0 * @settings[:numerical_hessian_step]**2) *
							(y[0] + y[1] - y[2] - y[3])
					when 4
						d2dxy = 1.0/(600.0 * @settings[:numerical_hessian_step]**2) *
							(-63.0 * (y[0] + y[1] + y[2] + y[3]) +
							  63.0 * (y[4] + y[5] + y[6] + y[7]) +
							  44.0 * (y[8] + y[9] - y[10] - y[11]) +
							  74.0 * (y[12] + y[13] - y[14] - y[15]))
					end

					results.hessian[i,j] = d2dxy
					results.hessian[j,i] = d2dxy
				}
			}
		end

		return results
	end

	def cleanup_e_only
		@energy_calc.cleanup

		if @what.include?(:gradient)
			@displaced_calcs_g.each{|array_d|
				array_d.each{|calc|
					calc.cleanup
				}
			}
		end

		if @what.include?(:hessian)
			@displaced_calcs.each{|array_j|
				array_j.each{|array_d|
					array_d.each{|calc|
						calc.cleanup
					}
				}
			}
		end
	end
end
