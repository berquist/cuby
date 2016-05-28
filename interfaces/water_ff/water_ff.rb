################################################################################
#
# Interface WaterFF
#
# Author: Jan Rezac
# Date created: 2012-04-11
# License: Cuby4 license
# Description: Simple forcefield for water
# Status: Works
#
################################################################################

#===============================================================================
# Water forcefield: 3 center, flexible SPC/Fw model
# Journal of Chemical Physics 124, 024503 (2006)
#===============================================================================

module InterfaceWaterFf
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient]
	MODIFIER = false
	#=======================================================================


	def prepare_interface
		@charges_oh = {:O => -0.82, :H => 0.41}
		@sigma_o = 3.165492
		@epsilon_o = 0.1554253
		@r_oh = 1.012
		@a_hoh = 113.24 / 180.0 * Math::PI
		@k_oh = 1059.162
		@k_hoh = 75.90

		# Check geometry
		Cuby::error("At least 3 atoms are needed for calculation of water", self) if @geometry.size < 3
		Cuby::error("Number of atoms in water molecule or cluster must be divisible by 3", self) if @geometry.size % 3 != 0

		(@geometry.size / 3).times{|m|
			m = m * 3
			# Check atom order
			unless @geometry[m].element == :O && @geometry[m+1].element == :H && @geometry[m+2].element == :H
				Cuby::error("Water molecules in the geometry must be ordered in OHHOHH... manner", self)
			end
			# Check distances
			Cuby::error("OH distance too large (atoms #{m+1} and #{m+2})", self) if @geometry[m].distance(@geometry[m+1]) > 1.5
			Cuby::error("OH distance too large (atoms #{m+1} and #{m+3})", self) if @geometry[m].distance(@geometry[m+2]) > 1.5
		}
	end

	def calculate_interface
		return water_ff_energy_grad(@what.include?(:gradient))
	end

	def water_ff_energy_grad(grad = false)
		results = Results.new

		if grad
			results.gradient = Gradient.zero(@geometry.size)
		end

		# Bonds and angles
		energy_b = 0.0
		energy_a = 0.0
		(@geometry.size / 3).times{|m|
			m = m*3
			# bond 1
			dist = @geometry[m].distance(@geometry[m+1])
			energy_b += 0.5 * @k_oh * (dist - @r_oh)**2
			if grad
				deriv = @k_oh * (dist - @r_oh)
				results.gradient.add_internal_dist(@geometry, deriv, m, m+1)
			end

			# bond 2
			dist = @geometry[m].distance(@geometry[m+2])
			energy_b += 0.5 * @k_oh * (dist - @r_oh)**2
			if grad
				deriv = @k_oh * (dist - @r_oh)
				results.gradient.add_internal_dist(@geometry, deriv, m, m+2)
			end

			# angle H-O-H
			angle = Coordinate.angle(@geometry[m+1], @geometry[m], @geometry[m+2])
			energy_a += 0.5 * @k_hoh * (angle - @a_hoh)**2
			if grad
				deriv = @k_hoh * (angle - @a_hoh)
				results.gradient.add_internal_angle(@geometry, deriv, m+1, m, m+2)
			end
		}
		
		# Electrostatics
		energy_e = 0.0
		@geometry.each_index{|i|
			i.times{|j|
				next if i/3 == j/3 # No interaction within one molecule
				dist = @geometry[i].distance(@geometry[j])
				q1 = @charges_oh[@geometry[i].element]
				q2 = @charges_oh[@geometry[j].element]
				energy_e += q1 * q2 / dist / ANGSTROM2BOHR * HARTREE2KCAL
				if grad
					deriv = -1.0 * q1 * q2 / dist**2 / ANGSTROM2BOHR * HARTREE2KCAL
					results.gradient.add_internal_dist(@geometry, deriv, i, j)
				end
			}
		}

		# L-J potential
		energy_lj = 0.0
		(@geometry.size / 3).times{|m|
			m.times{|n|
				r = @geometry[m*3].distance(@geometry[n*3])
				x6 = (@sigma_o / r)**6
				energy_lj += 4.0 * @epsilon_o * (x6**2 - x6)
				if grad
					deriv = 4.0 * @epsilon_o * (-12.0 * @sigma_o**12 / r**13  +  6.0 * @sigma_o**6 / r**7)
					results.gradient.add_internal_dist(@geometry, deriv, m*3, n*3)
				end


			}
		}

		results.energy = energy_b + energy_a + energy_e + energy_lj

		return results
	end
end
