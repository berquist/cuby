################################################################################
#
# Waterball interface
#
# Author: Jan Rezac
# Date created: 2014-05-27
# License: Cuby4 license
# Description: Potential holding a water droplet around a solute
# Status: Works
#
################################################################################

# Potential of Beglov & Roux[1]: Uses adaptive sphere radius to maintain
# a constant pressure of 1 atmosphere. The implementation here is only partial,
# the electrostatic interaction of the solvent with the cavity is not
# calculated, nor is the angular correction described in the paper. Even this
# partial implementation yields results sufficient for use, and is pricipaly
# similar to the hard sphere model of Karplus[2]
# [1] D. Beglov, B. Roux, Finite representation of an infinite bulk system: solvent boundary potential for computer simulations, Journal of Chemical Physics 100 (12), 9050-9063, 1994
# [2] A. Brunger, C. L. Brookds III, M. Karplus, Stochastic boundary condition for molecular dynamics simulations of ST2 water, Chemical Physics Letters 105(5), 495-498, 1984

### To be done:
### + Add angular correction
### + Add empirical correction for density
### + Add electrostatic term


module InterfaceWaterBall
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "Implemented parts work well"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient]
	MODIFIER = true
	#=======================================================================


	def prepare_interface
		# Ball center: :origin or :solute center (selection)
		@waterball_center = @settings[:waterball_center] # :origin or :solute
		if @waterball_center == :solute || @waterball_center == :fixed_solute
			@waterball_solute_list = @geometry.atomlist_from_selection(@settings[:waterball_solute])
			@waterball_solute_geo = @geometry.geometry_from_list(@waterball_solute_list)

			# Check
			if @waterball_solute_list.size < 10
				Cuby::recommendation("The solute defining the waterball is small, consider using larger one;
the force acting on the solute can be too large (it is distributed among all
the atoms of the solute).")
			end
		end

		# Make list of all water oxygens (atomlist)
		@waterball_oxygens = @geometry.atomlist_from_selection(":WAT@O")
		Cuby::log.puts_debug "#{@waterball_oxygens.size} water oxygens found\n\n"

		# Parameters
		@waterball_dvdw = 2.6
	end

	def calculate_interface
		return waterball_calculate(what.include?(:gradient))
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def waterball_calculate(grad)
		# Intialize results
		results = Results.new
		results.energy = 0.0
		results.gradient = Gradient.zero(@geometry.size) if grad

		# Get current center
		center = waterball_current_center
		Cuby::log.puts_debug "Center: #{center}"

		# Get Rmax and Imax
		r_max = 0.0
		i_max = nil
		@waterball_oxygens.each {|i|
			r = center.distance(@geometry[i])
			if r > r_max
				r_max = r 
				i_max = i
			end
		}

		# Parameters - from original paper
		delta_r_vdw = @waterball_dvdw
		x_c = 15.393
		a0 = -1.665
		a1 = 0.562
		a2 = -0.0728
		a3 = 0.00426
		a4 = -0.0000925
		a_c = 0.084
		b1 = 1.320
		b2 = -0.841
		b3 = -0.00160
		b4 = -8.37529889181283
		b_c = -8.475
		b5 = -9.56451953390558
		b6 = 1.6

		r_measure = 6.0
		n_measure = 0
		count = 0
		rmaxgrad = 0.0
		@waterball_oxygens.each {|i|
			r = center.distance(@geometry[i])
			r_vdw = r_max + delta_r_vdw
			x = r - r_vdw

			if grad
				direction = (@geometry[i] - center).normalize
				gradsize = 0.0
			end

			if r_vdw < x_c
				results.energy += a0 + a1 * r_vdw + a2 * r_vdw**2 + a3 * r_vdw**3 + a4 * r_vdw**4
				if grad
					de_drmax = 4.0*a4*r_vdw**3 + 3.0*a3*r_vdw**2 + 2.0*a2*r_vdw + a1
					rmaxgrad += de_drmax
				end
			else
				results.energy += a_c
				# no gradient
			end

			if x <= -5.0
				results.energy += b_c
				# no gradient
			elsif x <= 0.0
				results.energy += 1.0/(b2*(1.0+x**2/b1)) + b3 * x**2 + b4
				if grad
					de_dr = 2.0*b3*x - (2.0*x)/(b1*b2*(x**2/b1+1.0)**2)
					gradsize += de_dr
					rmaxgrad += -de_dr
				end
			else
				results.energy += b5 + b6 * x**2
				if grad
					de_dr = 2.0*b6*x
					gradsize += de_dr
					rmaxgrad += -de_dr
				end
			end

			if grad
				results.gradient[i] = direction * gradsize if gradsize != 0.0
			end

			n_measure += 1 if r < r_measure # count atoms for density measurement

			# Count shell molecules
			if r > r_max - 2.0 # outer shell molecules
				count += 1
			end
		}

		# Gradient on the atom defining r_max
		if grad
			direction = (@geometry[i_max] - center).normalize
			results.gradient[i_max] += direction * rmaxgrad
		end

		Cuby::log.puts_debug "Rmax: #{r_max}"
		Cuby::log.puts_debug "Molecules in outer shell: #{count}"

		# Calculate density
		volume = 4.0/3.0 * Math::PI * r_measure**3
		Cuby::log.puts_debug "Density #{n_measure / volume} A^-3"

		# Energy from pressure and surface tension
		surface = 4.0 * Math::PI * r_max**2
		volume = 4.0/3.0 * Math::PI * r_max**3
		# Pa = 6947694303.36359  * program units
		# 1 atm = 1.4584007e-05
		a = 1.4584007e-05 * volume # pV at 1 atm
		b = 0.033 * surface # surf. ten * surface
		results.energy += a + b

		# Gradient from pressure and surface tension
		if grad
			grad = (3.0 * a + 2.0 * b) / r_max
			results.gradient[i_max] += grad * (@geometry[i_max] - center).normalize
		end

		if @waterball_center == :fixed_solute
			waterball_solute_potential(results, grad)
		elsif @waterball_center == :solute && grad
			# Opposite gradient should be distributed onto the solute
			avggrad = Coordinate[0.0,0.0,0.0]
			results.gradient.each{|g|
				avggrad += g
			}
			avgggrad = avggrad / @waterball_solute_list.size
			@waterball_solute_list.each{|i|
				results.gradient[i] -= avgggrad
			}
		end

		return results
	end
	
	def waterball_current_center # => Coordinate
		#: Returns center of the sphere according to current setup
		case @waterball_center
		when :origin
			return Coordinate[0.0,0.0,0.0]
		when :fixed_solute
			return Coordinate[0.0,0.0,0.0]
		when :solute
			return @waterball_solute_geo.center
		end
	end

	def waterball_solute_potential(results, grad)
		#: Harmonic potential applied to solute when :fixed_solute used
		k = @settings[:waterball_solute_k]
		solute_center = @waterball_solute_geo.center
		center = waterball_current_center
		r = center.distance(solute_center)
		results.energy += 0.5 * k * r**2
		if grad
			g = k * r
			direction = (solute_center - center).normalize
			@waterball_solute_list.each { |i|
				results.gradient[i] += g * direction * (1.0 / @waterball_solute_list.size)
			}
		end
	end
end
