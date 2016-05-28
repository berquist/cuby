################################################################################
#
# Module Func
#
# Status: Works
#
################################################################################

#: This module serves as a conainer for several clases representing the functions.
#: Instance of the class stores all the parameters, only the most important variable
#: is provided when the value is calculated.


module Func

	class PolySwitching
		# Polynomial switching function with smooth first and second derivative

		def initialize (r_0, r_1, direction = :sw_0_to_1)
			@r0 = r_0
			@r1 = r_1
			@direction = direction
		end 

		def calc(r)
			if @direction == :sw_0_to_1
				return 0.0 if r <= @r0
				return 1.0 if r >= @r1
				interval = @r1 - @r0
				x = (r - @r0) / interval
				return -20.0*x**7 + 70.0*x**6 - 84.0*x**5 + 35.0*x**4
			else
				return 1.0 if r <= @r0
				return 0.0 if r >= @r1
				interval = @r1 - @r0
				x = (r - @r0) / interval
				return 1.0 - (-20.0*x**7 + 70.0*x**6 - 84.0*x**5 + 35.0*x**4)
			end
		end

		def deriv(r)
			return 0.0 if r <= @r0
			return 0.0 if r >= @r1
			interval = @r1 - @r0
			x = (r - @r0) / interval
			if @direction == :sw_0_to_1
				return (-140.0*x**6 + 420.0*x**5 - 420.0*x**4 + 140.0*x**3) / interval
			else
				return -(-140.0*x**6 + 420.0*x**5 - 420.0*x**4 + 140.0*x**3) / interval
			end

		end

		def deriv2(r)
			return 0.0 if r <= @r0
			return 0.0 if r >= @r1
			interval = @r1 - @r0
			x = (r - @r0) / interval
			if @direction == :sw_0_to_1
				return (-840.0*x**5 + 2100.0*x**4 - 1680.0*x**3 + 420.0*x**2) / interval
			else
				return -(-840.0*x**5 + 2100.0*x**4 - 1680.0*x**3 + 420.0*x**2) / interval
			end
		end
	end

	class DampingFunction
		#: Smooth damping function. Type selects the form:
		#* :df_1_to_0 => ^^^\___ (default)
		#* :df_0_to_1 => ___/^^^

		attr_accessor :r0, :exponent, :type

		def initialize(r0,exponent,type = :df_1_to_0)
			#: Creates new instance of DampingFunction
			@r0 = r0
			@exponent = exponent
			@type = type
		end

		def calc(r)
			#: Calculates value of the damping function at point r
			case @type
			when :df_0_to_1
				# 0 __/^^ 1
				return 1.0 / (1.0 + Math::exp(-@exponent * (r/@r0 - 1.0)))
			when :df_1_to_0
				# 1 ^^\__ 0
				return 1.0 - 1.0 / (1.0 + Math::exp(-@exponent * (r/@r0 - 1.0)))
			end
		end

		def deriv(r)
			d = 1.0 / (1.0 + Math::exp(-@exponent*(r/@r0-1.0)))**2 * @exponent/@r0 * Math::exp(-@exponent*(r/@r0-1.0))
			case @type
			when :df_0_to_1
				return d
			when :df_1_to_0
				return -1.0 * d
			end
		end

		def deriv2(r)
			d2 = (2.0*@exponent**2*Math::exp(-2.0*@exponent*(r/@r0-1.0))) / 
				(@r0**2*(Math::exp(-@exponent*(r/@r0-1.0))+1.0)**3) - 
				(@exponent**2*Math::exp(-@exponent*(r/@r0-1.0))) /
				(@r0**2*(Math::exp(-@exponent*(r/@r0-1.0))+1.0)**2)
			case @type
			when :df_0_to_1
				return d2
			when :df_1_to_0
				return -1.0 * d2
			end
		end
	end
end
