################################################################################
#
# Class Thermodynamics
#
# Author: Jan Rezac
# Date created: 2009-07-16
# License: Cuby license
# Description: Thermodynamic equations
# Status: Planned
#
################################################################################

#: Calculates thermodynamic properties using ideal gas approximation

class Thermodynamics

	attr_reader :pressure
	attr_reader :temperature
	attr_reader :symmetry_number

	#=======================================================================
	# Initialize
	#=======================================================================

	def initialize(settings, geometry, frequencies, freq_limit = nil)
		@geometry = geometry
		@frequencies = frequencies
		@freq_limit = freq_limit

		@temperature = settings[:thermo_temperature]
		@pressure = settings[:thermo_pressure]
		@symmetry_number = settings[:thermo_symmetry_number]
		@thermo_low_mode_fix = settings[:thermo_low_mode_fix]


		# Precalculate average moment of inertia
		if @geometry.size < 3
			@inertia_av = nil
		else
			moments = @geometry.moment_of_inertia_diagonalized[1].to_a # in mass_unit * angstrom^2
			# Moment of inertia in g/mol * angstrom^2
			@moments_of_inertia_gmola =  (Vector.from_array(moments) * MASSUNIT2KG * 1000 * AVOGADRO).to_a
			@inertia_av = 0.0
			moments.each{|i| @inertia_av += i * MASSUNIT2KG * ANGSTROM2M**2 }
			@inertia_av = @inertia_av / 3.0
		end
	end

	#=======================================================================
	# Setters
	#=======================================================================
	# Using separate setters could be used to cache some partial results
	
	def temperature=(temperature)
		@temperature = temperature
	end

	def pressure=(pressure)
		@pressure = pressure
	end

	def symmetry_number=(symmetry_number)
		@symmetry_number = symmetry_number
	end

	
	#=======================================================================
	# Partition function
	#=======================================================================
	
	def q_translation
		# Ideal gas equation:
		volume = BOLTZMANN_CONSTANT_SI*@temperature/@pressure
		# Mass (kg/mol)
		mass = @geometry.mass * MASSUNIT2KG
		return (2.0 * Math::PI * mass * BOLTZMANN_CONSTANT_SI * @temperature / PLANCK_CONSTANT_SI**2)**1.5 * volume
	end

	def dlnq_dT_translation
		return 1.5/@temperature
	end

	def q_rotation
		case @geometry.size
		when 1
			# Atom
			return 1.0
		when 2
			# Linear molecule
			com = @geometry.center_of_mass
			moment = @geometry[0].distance(com)**2 * @geometry[0].mass +
				@geometry[1].distance(com)**2 * @geometry[1].mass
			moment *= MASSUNIT2KG * ANGSTROM2M**2
			rottemp = PLANCK_CONSTANT_SI**2 / (8.0 * Math::PI**2 * moment * BOLTZMANN_CONSTANT_SI)
			return 1.0 / @symmetry_number * @temperature / rottemp
		else
			# Polyatomic molecule
			moments = @geometry.moment_of_inertia_diagonalized[1].to_a # in mass_unit * angstrom^2
			moments = moments.map{|i| i * MASSUNIT2KG * ANGSTROM2M**2 }
			rottemps = moments.map{|i|
				PLANCK_CONSTANT_SI**2 / (8.0 * Math::PI**2 * i * BOLTZMANN_CONSTANT_SI)
			}
			return Math::PI**0.5 / @symmetry_number *@temperature**1.5 /
				(rottemps[0] * rottemps[1] * rottemps[2])**0.5
		end
	end

	def dlnq_dT_rotation
		case @geometry.size
		when 1
			# Atom
			return 0.0
		when 2
			# Linear molecule
			return 1.0/@temperature
		else
			# Polyatomic molecule
			return 1.5/@temperature
		end
	end

	#=======================================================================
	# Final equations : entropy
	#=======================================================================
	
	def entropy_translation
		return GAS_CONSTANT_SI * (
			1.0 +
			Math::log(q_translation) +
			@temperature * dlnq_dT_translation
		)
	end

	def entropy_vibration_mode(freq_cm)
		freq = freq_cm * CM2HZ
		hfkt = freq * PLANCK_CONSTANT_SI / BOLTZMANN_CONSTANT_SI / @temperature
		
		return GAS_CONSTANT_SI * (
			hfkt / (Math.exp(hfkt) - 1.0) -
			Math.log(1.0-Math.exp(-hfkt))
		)
	end

	def entropy_vibration_mode_low_fix(freq_cm)
		freq = freq_cm * CM2HZ
		hfkt = freq * PLANCK_CONSTANT_SI / BOLTZMANN_CONSTANT_SI / @temperature
		
		# Original vibrational entropy
		vibrational_s = GAS_CONSTANT_SI * (
			hfkt / (Math.exp(hfkt) - 1.0) -
			Math.log(1.0-Math.exp(-hfkt))
		)


		# Moment of inertia of a rotor with the same frequency
		mu = PLANCK_CONSTANT_SI / 8.0 / Math::PI**2 / freq

		# Average molecular moment of inertia
		if @inertia_av
			b_av = @inertia_av
		else
			Cuby::error("Low mode vibration correction can not be used form molecules with less than 3 atoms")
		end

		# Reduced moment of inertia
		mu_r = mu*b_av/(mu + b_av)

		# 'Rotational' entropy
		rotational_s = GAS_CONSTANT_SI * (
			0.5 + Math.log(
				(8.0 * Math::PI**3 * mu_r * BOLTZMANN_CONSTANT_SI * @temperature / PLANCK_CONSTANT_SI**2)**0.5
			)
		)

		# Weighting function
		f0 = 100 * CM2HZ
		w = 1.0 / (1.0 + (f0 / freq)**4)

		# Weighted mix of entropies
		s = w * vibrational_s + (1.0-w) * rotational_s

		return s
	end

	def entropy_vibration
		s = 0.0
		@frequencies.each_index{|i|
			f = @frequencies[i]
			return s if i == @freq_limit
			s += entropy_vibration_mode(f) if f > 0
		}
		return s
	end

	def entropy_vibration_low_fix(invert_negative = false)
		s = 0.0
		@frequencies.each_index{|i|
			f = @frequencies[i]
			return s if i == @freq_limit
			if f > 0
				s += entropy_vibration_mode_low_fix(f) 
			elsif f < 0 && invert_negative
				s += entropy_vibration_mode_low_fix(-f)
			end
		}
		return s
	end

	def entropy_rotation
		return GAS_CONSTANT_SI * (
			Math::log(q_rotation) +
			@temperature * dlnq_dT_rotation
		)
	end

	#=======================================================================
	# Final equations : energy
	#=======================================================================
	
	def energy_translation
		return GAS_CONSTANT_SI * @temperature * 1.5
	end
	
	def energy_rotation
		return GAS_CONSTANT_SI * @temperature**2 * dlnq_dT_rotation
	end

	def energy_vibration_mode(freq_cm)
		freq = freq_cm * CM2HZ
		hfk = freq * PLANCK_CONSTANT_SI / BOLTZMANN_CONSTANT_SI
		return GAS_CONSTANT_SI * hfk * (0.5 + 1.0/(Math.exp(hfk/@temperature) - 1.0))
	end

	def energy_vibration
		e = 0.0
		@frequencies.each_index{|i|
			f = @frequencies[i]
			return e if i == @freq_limit
			e += energy_vibration_mode(f) if f > 0
		}
		return e
	end

	#=======================================================================
	# Final equations : Cp
	#=======================================================================
	
	#=======================================================================
	# Final equations : Cv
	#=======================================================================

	#=======================================================================
	# Printing
	#=======================================================================
	
	def print(energy)
		Cuby::log.puts "Temperature: #{temperature} K"
		Cuby::log.puts
		if @moments_of_inertia_gmola
			Cuby::log.puts "Moments of inertia:"
			Cuby::log.puts "   Ia #{'%12.3f' % (@moments_of_inertia_gmola[0])} g/mol*Angstrom"
			Cuby::log.puts "   Ib #{'%12.3f' % (@moments_of_inertia_gmola[1])} g/mol*Angstrom"
			Cuby::log.puts "   Ic #{'%12.3f' % (@moments_of_inertia_gmola[2])} g/mol*Angstrom"
			Cuby::log.puts
		end
		s_t = entropy_translation * J2CAL
		s_r = entropy_rotation * J2CAL
		case @thermo_low_mode_fix
		when :none
			# original
			s_v = entropy_vibration * J2CAL
		when :as_rotations
			# fix
			s_v = entropy_vibration_low_fix * J2CAL
		when :as_rotations_including_negative
			# fix
			s_v = entropy_vibration_low_fix(true) * J2CAL
		end
		Cuby::log.puts "S(trans): #{'%8.3f' % (s_t)} cal/mol/K"
		Cuby::log.puts "S(rot):   #{'%8.3f' % (s_r)} cal/mol/K"
		Cuby::log.puts "S(vib):   #{'%8.3f' % (s_v)} cal/mol/K"
		Cuby::log.puts "sum S:    #{'%8.3f' % (s_v + s_r + s_t)} cal/mol/K"
		Cuby::log.puts

		e_t = energy_translation * J2KCAL
		e_r = energy_rotation * J2KCAL
		e_v = energy_vibration * J2KCAL
		Cuby::log.puts "E(trans): #{'%8.3f' % (e_t)} kcal/mol"
		Cuby::log.puts "E(rot):   #{'%8.3f' % (e_r)} kcal/mol"
		Cuby::log.puts "E(vib):   #{'%8.3f' % (e_v)} kcal/mol (includes ZPVE)"
		Cuby::log.puts "sum E:    #{'%8.3f' % (e_t + e_r + e_v)} kcal/mol"
		Cuby::log.puts
		ts = (s_v + s_r + s_t) * temperature / 1000
		e_temp = e_t + e_r + e_v + energy
		Cuby::log.puts "E(T) = E(el) + E(t,r,v)"
		Cuby::log.puts "G = E(T) - TS"
		Cuby::log.puts "E(el) =   #{'%8.3f' % (energy)} kcal/mol"
		Cuby::log.puts "E(T)  =   #{'%8.3f' % (e_temp)} kcal/mol"
		Cuby::log.puts "TS    =   #{'%8.3f' % (ts)} kcal/mol"
		Cuby::log.puts "G     =   #{'%8.3f' % (e_temp - ts)} kcal/mol"
	end
end
