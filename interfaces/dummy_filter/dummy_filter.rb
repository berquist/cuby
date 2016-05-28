################################################################################
#
# Dummy atom filter interface
#
# Author: Jan Rezac
# Date created: 2014-09-25
# License: Cuby4 license
# Description: Filter for dummy atoms
# Status: Works
#
################################################################################

#===============================================================================
# Removes dummy atoms from the geometry and adds them back to results
#===============================================================================

module InterfaceDummyFilter
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
		InputBlock[:calculation, :required, "Setup for the calculation"]
	]
	#=======================================================================

	def prepare_interface
		# Child settings
		calc_settings = @settings.block(:calculation)

		# New geometry without dummy atoms
		geometry = filtered_geometry(@geometry)

		# Build calculation
		@calc = Calculation.new(@name+'_filtered', calc_settings, geometry)

		# Prepare calculation
		@calc.prepare(@what)

	end

	def queue_interface(queue)
		@calc.send_to_queue(queue)
	end

	def compose_interface
		results = Results.new

		results.energy = @calc.results.energy
		results.energy_components = @calc.results.energy_components
		
		if @what.include?(:gradient)
			results.gradient = Gradient.new
			i_real = 0
			@geometry.each_index{|i|
				atom = @geometry[i]
				if atom.element == :X
					results.gradient << Coordinate[0,0,0]
				else
					results.gradient << @calc.results.gradient[i_real]
					i_real += 1
				end
			}
		end

		return results
	end

	def cleanup_interface
		@calc.cleanup
	end

	def filtered_geometry(geo)
		ng = Geometry.new
		geo.each{|atom|
			ng << atom unless atom.element == :X
		}
		return ng
	end
end
