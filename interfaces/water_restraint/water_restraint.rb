################################################################################
#
# Water restraint interface
#
# Author: Jan Rezac
# Date created: 2014-05-23
# License: Cuby4 license
# Description: Restraint preventing mixing QM and MM water in QM/MM MD
# Status: Works
#
################################################################################


module InterfaceWaterRestraint
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "Passed tests, gradient correct"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient]
	MODIFIER = true
	#=======================================================================


	def prepare_interface
		# Get sphere center
		@sphere_center_list = @geometry.atomlist_from_selection(@settings[:qmmm_wat_center])
		@sphere_center_geo = @geometry.geometry_from_list(@sphere_center_list)
		@sphere_n = @settings[:qmmm_wat_num]
		@k = @settings[:qmmm_wat_k]

		# Make a list of water molecules
		@connectivity = Connectivity.new(geometry)
		@water_oxygens = find_water_oxygens(@geometry, @connectivity)
		Cuby::log.puts_debug("Found #{@water_oxygens.count} water molecules")
		# Sort oxygens by distance from center
		sort_waters(@water_oxygens, @sphere_center_geo.center)

		# Select N closest waters
		@qmwaters = @water_oxygens[0..@sphere_n-1]
		
		# For debugging, mix the waters: 
		# @qmwaters = [0,1,2,3, 4,5,8,10].map{|i| @water_oxygens[i]}
	end

	def calculate_interface
		return wr_calculate(what.include?(:gradient))
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def wr_energy(sphere_com, active_mmwaters, active_qmwaters, geometry, silent = false)
		# find closest MM water and most distant QM water
		firstmm = active_mmwaters.first
		lastqm = active_qmwaters.last

		r_mm = (sphere_com - geometry[firstmm]).abs
		r_qm = (sphere_com - geometry[lastqm]).abs


		# Energy: Harmonic potential
		energy = 0.0
		active_mmwaters.each{|i|
			ri = (sphere_com - geometry[i]).abs
			energy += 0.5 * @k * (ri - r_qm)**2

		}
		active_qmwaters.each{|i|
			ri = (sphere_com - geometry[i]).abs
			energy += 0.5 * @k * (ri - r_mm)**2

		}

		Cuby.log.puts_debug("Water restraint applied, r_mm = #{'%.3f' % r_mm}, r_qm = #{'%.3f' % r_qm}, E = #{'%.3f' % energy}") unless silent
		return energy
	end
	
	def wr_calculate(grad)
		# Intialize results
		results = Results.new
		results.energy = 0.0
		results.gradient = Gradient.zero(@geometry.size) if grad

		# Sort waters, pick the closest ones
		sphere_com = @sphere_center_geo.center
		sort_waters(@water_oxygens, sphere_com)
		sort_waters(@qmwaters, sphere_com)
		closewaters = @water_oxygens[0..@sphere_n-1]

		# Nothing to do when all QM waters are the closest ones
		if (@qmwaters - closewaters).size == 0
			Cuby.log.puts_debug("Water restraint not active")
			return results
		end

		# Get active waters
		active_mmwaters = []
		active_qmwaters = []
		count = 0
		i = 0
		mm_found = false
		while count < @sphere_n do
			if @qmwaters.include?(@water_oxygens[i])
				# QM water
				active_qmwaters << @water_oxygens[i] if mm_found
				count += 1
			else
				# MM water
				mm_found = true
				active_mmwaters << @water_oxygens[i]
			end
			i += 1
		end

		# Energy		
		results.energy = wr_energy(sphere_com, active_mmwaters, active_qmwaters, @geometry)

		# Gradient applies to more atoms, the QM and MM water and the solute
		if grad
			displacement = 0.0001
			gradatoms = @sphere_center_list + active_mmwaters + active_qmwaters
			g2 = @geometry.deep_copy # Calulate on a duplicate geometry, to not to disturb the original one
			sphere_center_geo2 = g2.geometry_from_list(@sphere_center_list)
			gradatoms.each{|i|
				atom = g2[i]
				fgrad = Coordinate.new
				3.times{|c|
					# Plus d
					atom[c] += displacement
					# For sphere center atoms, com must be recalculated
					sphere_com2 = sphere_center_geo2.center
					ep = wr_energy(sphere_com2, active_mmwaters, active_qmwaters, g2, true)
					# Minus d
					atom[c] -= displacement * 2
					# For sphere center atoms, com must be recalculated
					sphere_com2 = sphere_center_geo2.center
					em = wr_energy(sphere_com2, active_mmwaters, active_qmwaters, g2, true)
					# Restore coordinates
					atom[c] += displacement
					# Gradient
					fgrad[c] = (ep - em) / (displacement * 2.0)
				}
				results.gradient[i] += fgrad
			}

		end

		return results
	end

	def find_water_oxygens(geometry, connectivity) # => Array of atom indices
		oxygens = []
		geometry.each_with_index{|atom,i|
			if atom.element == :O &&
			   connectivity.hybridization(i) == 2 &&
			   geometry[connectivity.bound_atoms(i)[0]].element == :H &&
			   geometry[connectivity.bound_atoms(i)[1]].element == :H &&
				# Water oxygen, save atom index
				oxygens << i
			end
		}
		return oxygens
	end

	def sort_waters(water_oxygens, sphere_com)
		water_oxygens.sort!{|ia, ib| @geometry[ia].distance(sphere_com) <=> @geometry[ib].distance(sphere_com)}
	end

end
