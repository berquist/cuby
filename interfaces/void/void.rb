################################################################################
#
# Void interface
#
# Author: Jan Rezac
# Date created: 2012-04-11
# License: Cuby4 license
# Description: Interface returning zero values
# Status: Works
#
################################################################################

#===============================================================================
# Interface returning zero values of energy, gradient and hessian.
# Useful for testing.
#===============================================================================

module InterfaceVoid
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient, :hessian]
	MODIFIER = true
	#=======================================================================


	def prepare_interface
		# Nothing to be done
	end

	def calculate_interface
		results = Results.new

		if @what.include?(:energy)
			results.energy = 0.0
		end

		if @what.include?(:gradient)
			results.gradient = Gradient.zero(@geometry.size)
		end

		if @what.include?(:hessian)
			results.hessian = Hessian.zero(@geometry.size * 3)
		end

		sleep(@settings[:void_delay])

		return results
	end
end
