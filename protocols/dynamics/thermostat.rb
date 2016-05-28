module Thermostat
	#=======================================================================
	# No thermostat, plain NVE
	#=======================================================================
	module None
		attr_reader :thermostat_temp

		def thermostat_setup
		end

		def thermostat_set_t(temp)
		end

		def integrator
			integrator_verlet
		end

		def thermostat(t)
		end
	end

	#=======================================================================
	# Berendsen thermostat
	#=======================================================================
	module Berendsen
		attr_reader :thermostat_temp

		def thermostat_setup
			@thermostat_temp = @settings[:temperature]
			@berendsen_tau_t = @settings[:thermostat_tc]
		end

		def thermostat_set_t(temp)
			@thermostat_temp = temp
		end

		def integrator
			integrator_verlet
		end

		def thermostat(t = temperature)
			factor = Math.sqrt(1 + (@thermostat_temp/t-1.0) * @dt/@berendsen_tau_t)
			if factor < 0.8
				factor = 0.8
			end
			if factor > 1.25
				factor = 1.25
			end
			rescale_velocities(factor)
		end
	end

	#=======================================================================
	# Andersen thermostat
	#=======================================================================
	module Andersen
		attr_reader :thermostat_temp

		def thermostat_setup
			tc = @settings[:thermostat_tc]
			@thermostat_temp = @settings[:temperature]
			@andersen_tc = @dt / @geometry.size / tc
		end

		def thermostat_set_t(temp)
			@thermostat_temp = temp
		end

		def integrator
			integrator_verlet
		end

		def thermostat(t = nil)
			@geometry.each_index {|atom_i|
				if rand < @andersen_tc
					generate_maxwell_velocity(atom_i, @thermostat_temp)
				end
			}
		end
	end

	#=======================================================================
	# Nose-Hoover thermostat
	#=======================================================================
	#: Formulation from Gromacs manual
	
	module NoseHoover
		attr_reader :thermostat_temp

		def thermostat_setup
			# setup
			@thermostat_temp = @settings[:temperature]
			# init
			@nose_hoover_xi = 0.0
			@nose_hoover_vxi = 0.0
			@nose_hoover_q = @settings[:thermostat_tc] * @thermostat_temp / 4.0 / Math::PI**2
		end

		def thermostat_set_t(temp)
			@thermostat_temp = temp
			@nose_hoover_q = @settings[:thermostat_tc] * @thermostat_temp / 4.0 / Math::PI**2
		end

		def integrator
			# Using special integrator here...
			selection_grad = selection_gradient
			@geometry.each_with_index{|atom, i|
				tempcoord = atom.to_coordinate
				atom.plus!(atom - @oldcoord[i] - selection_grad[i] * (@dt2 / atom.mass))
				atom.minus!((atom - @oldcoord[i]) * @dt * @nose_hoover_xi) # This line is added to normal Verlet
				@oldcoord[i] = tempcoord
			}
			return nil
		end

		def thermostat(t = temperature)
			@nose_hoover_vxi = 1.0 / @nose_hoover_q * (t - @thermostat_temp)
			# propagate xi
			@nose_hoover_xi += @nose_hoover_vxi *  @dt
		end
	end

	#=======================================================================
	# Bussi thermostat
	#=======================================================================
	#: Bussi, Donadio, Parrinello, J. Chem. Phys. 126, 014101 (2007)
	
	#: Global rescaling factor applied to velocities (like Berendsen), with
	#: randomness added into that factor to ensure sampling of canonical
	#: ensemble (unlike Berendsen, and better).
	
	module Bussi
		attr_reader :thermostat_temp

		def thermostat_setup
			@thermostat_temp = @settings[:temperature]
			@bussi_tau_t = @settings[:thermostat_tc]
		end

		def thermostat_set_t(temp)
			@thermostat_temp = temp
		end

		def integrator
			integrator_verlet
		end

		def thermostat(t = temperature)
			# get kinetic energies
			k_set = @thermostat_temp / (2.0 / (3.0 * @geometry.size * BOLTZMANN))
			k = t / (2.0 / (3.0 * @geometry.size * BOLTZMANN))

			nf = @geometry.size * 3.0
			
			e = Math.exp(-@dt/@bussi_tau_t)
			r1 = rand_gauss
			sum_gauss_sq = 0.0
			(nf-1).to_i.times {sum_gauss_sq += rand_gauss**2}

			a2 = 	e + 
				k_set / (k * nf) * (1.0 - e) * (r1**2 + sum_gauss_sq) +
				2.0 * Math.exp(-@dt/@bussi_tau_t/2.0) * Math.sqrt(k_set / (k * nf) * (1.0 - e)) * r1

			factor = Math.sqrt(a2)

			rescale_velocities(factor)
		end

		@@rand_gauss_saved = nil
		def rand_gauss
			# algorithm from http://www.taygeta.com/random/gaussian.html

			# return saved value if it exists
			if @@rand_gauss_saved
				retval = @@rand_gauss_saved
				@@rand_gauss_saved = nil
				return retval
			end

			# or calculate new pair
			w = 10.0
			x1 = x2 = 0.0
			while w >= 1.0
				x1 = 2.0 * rand - 1.0
				x2 = 2.0 * rand - 1.0
				w = x1 * x1 + x2 * x2
			end
			w = Math.sqrt( (-2.0 * Math.log( w ) ) / w )

			# save second result
			@@rand_gauss_saved = x2 * w

			return x1 * w 
		end
	end

	#=======================================================================
	# Langevin dynamics
	#=======================================================================
	#--------------------------------------------------
	# module Langevin
	# 	def thermostat_setup
	# 	end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	def integrator(t = temperature)
	# 		# Using special integrator here...
	# 	end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	def thermostat
	# 	end
	# end
	#-------------------------------------------------- 
end
