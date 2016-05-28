################################################################################
#
# Interface Electrostatics
#
# Author: Jan Rezac
# Date created: 2015-04-11
# License: Cuby4 license
# Description: Coulomb electostatics
# Status: Works
#
################################################################################

#===============================================================================
# Electrostatics based on atomic charges
#===============================================================================

module InterfaceElectrostatics
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
		# Load atomic charges from file
		if @settings.set?(:atomic_charges_read)
			if @settings[:atomic_charges_read].downcase == "from_geometry"
				@charges = AtomicCharges.from_geometry(@geometry)
			else
				@charges = AtomicCharges.from_file(@settings[:atomic_charges_read])
			end
		else
			Cuby::error "Atomic charges must be read from file specified by atomic_charges_read keyword"
		end

		if [:gaussian, :slater].include?(@settings[:electrostatics_model])
			# Initialize hubbard U table
			# Values from DFTB 3-OB parameter set
			@u = {
				:H => 0.4195,
				:C => 0.3647,
				:N => 0.4309,
				:O => 0.4954,
				:B => 0.2958,
				:S => 0.3288,
				:P => 0.2894,
				:F => 0.5584,
				:Cl => 0.3668,
				:Br => 0.327662,
				:I => 0.284213
			}

			@geometry.elements_in_system.each{|e|
				unless @u.has_key?(e)
					Cuby::error "Electrostatics model parameters not available for element #{e}"
				end
			}
		end
	end

	def calculate_interface
		case @settings[:electrostatics_model]
		when :point_charges
			return coulomb_energy_grad(@what.include?(:gradient))
		when :gaussian
			return gaussian_energy_grad(@what.include?(:gradient))
		when :slater
			return slater_energy_grad(@what.include?(:gradient))
		end
	end

	def coulomb_energy_grad(grad = false)
		results = Results.new

		if grad
			results.gradient = Gradient.zero(@geometry.size)
		end

		# Electrostatics
		energy_e = 0.0
		@geometry.each_index{|i|
			i.times{|j|
				dist = @geometry[i].distance(@geometry[j]) * ANGSTROM2BOHR
				q1 = @charges[i]
				q2 = @charges[j]
				energy_e += q1 * q2 / dist * HARTREE2KCAL
				if grad
					deriv = -1.0 * q1 * q2 / dist**2 / ANGSTROM2BOHR * HARTREE2KCAL
					results.gradient.add_internal_dist(@geometry, deriv, i, j)
				end
			}
		}

		results.energy = energy_e

		return results
	end

	def gaussian_energy_grad(grad = false)
		results = Results.new

		if grad
			raise "Gradient not implemented for gaussian densities"
			results.gradient = Gradient.zero(@geometry.size)
		end

		# Electrostatics
		energy_e = 0.0
		@geometry.each_index{|i|
			i.times{|j|
				dist = @geometry[i].distance(@geometry[j]) * ANGSTROM2BOHR
				q1 = @charges[i]
				q2 = @charges[j]
				fwhm_i = (8.0*Math::log(2)/Math::PI)**0.5 / @u[@geometry[i].element]
				fwhm_j = (8.0*Math::log(2)/Math::PI)**0.5 / @u[@geometry[j].element]
				c_ij = (4.0 * Math::log(2) / (fwhm_i**2 + fwhm_j**2))**0.5
				gamma = Math::erf(c_ij * dist) / dist
				energy_e += q1 * q2 * gamma * HARTREE2KCAL
				if grad
					deriv = 0.0
					results.gradient.add_internal_dist(@geometry, deriv, i, j)
				end
			}
		}

		results.energy = energy_e

		return results
	end

	def slater_energy_grad(grad = false)
		results = Results.new

		if grad
			raise "Gradient not implemented for slater densities"
			results.gradient = Gradient.zero(@geometry.size)
		end

		# Electrostatics
		energy_e = 0.0
		@geometry.each_index{|i|
			i.times{|j|
				dist = @geometry[i].distance(@geometry[j]) * ANGSTROM2BOHR
				q1 = @charges[i]
				q2 = @charges[j]
				ui = @u[@geometry[i].element]
				uj = @u[@geometry[j].element]
				energy_e += q1 * q2 * gamma(geometry[i].element, @geometry[j].element, ui, uj, dist) * HARTREE2KCAL
				if grad
					deriv = 0.0
					results.gradient.add_internal_dist(@geometry, deriv, i, j)
				end
			}
		}

		results.energy = energy_e

		return results
	end

	def gamma(element1, element2, u1, u2, r)
		a = 16.0/5.0*u1
		b = 16.0/5.0*u2

		# Same atom => average of Hubbard U's
		if r == 0
			if u1 == u2
				return 0.5 * (u1 + u2)
			else
				raise "Different Hubbard U on a single atom"
			end
		end

		# Calculate gamma using the S^g and S^f functions
		if u1 == u2
			sg = ((a**3 * r**3 + 9.0 * a**2 * r**2 + 33.0 * a * r + 48.0) * Math::exp(-a * r)) / (48.0 * r)
			g = 1.0 / r - sg
		else
			sf = ((a**4 * b) / (2.0 * (b**2 - a**2)**2) - (a**6 - 3.0 * a**4 * b**2) / 
			      ((b**2 - a**2)**3 * r)) * Math::exp(-b * r) +
			     ((a * b**4) / (2.0 * (a**2 - b**2)**2) - (b**6 - 3.0 * a**2 * b**4) / 
			      ((a**2 - b**2)**3 * r)) * Math::exp(-a * r)
			g = 1.0 / r - sf
		end

		return g
	end

end
