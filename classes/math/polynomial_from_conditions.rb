require "classes/algebra/algebra"
require "classes/math/polynomial"

# Constructs a polynomial satisfying a set of conditions.
# Order of the polynomial is equivalent to the number of conditions.

module PolynomialFromConditions
	include Algebra

	class Condition
		attr_accessor :type
		attr_accessor :x
		attr_accessor :value

		ALLOWED_TYPES = [:value, :first_deriv, :second_deriv, :third_deriv]

		def initialize(type, x, value = 0.0)
			unless ALLOWED_TYPES.include?(type)	
				raise "Unknown condition type"
			end

			@type = type
			@x = x
			@value = value
		end
	end

	def PolynomialFromConditions.build(conditions)
		#: Create polynomial that satisfies provided conditions.
		size = conditions.size
		m = Matrix.zero(size, size)
		my = Matrix.zero(size, 1)

		size.times{|i|
			x = conditions[i].x
			my[i,0] = conditions[i].value
			case conditions[i].type
			when :value
				 size.times{|j|
					 m[i,j] = x**(size - 1 - j)
				 }
			when :first_deriv
				 (size-1).times{|j|
					 m[i,j] = (size - 1 - j).to_f * x**(size - 2 - j)
				 }
			when :second_deriv
				 (size-2).times{|j|
					 m[i,j] = ((size - 1 - j)*(size - 2 - j)).to_f * x**(size - 3 - j)
				 }
			when :third_deriv
				 (size-3).times{|j|
					 m[i,j] = ((size - 1 - j)*(size - 2 - j)*(size - 3 - j)).to_f * x**(size - 4 - j)
				 }
			end
		}
		poly = (m.inverse * my).column_as_array(0)
		return  Polynomial.new(poly)
	end

end

