require "classes/algebra/algebra"

class Polynomial
	include Algebra

	def initialize(coefficients)
		#: Build polynomial from an array of coefficients (highest power first).
		@c = coefficients
	end

	def [](index)
		return @c[index]
	end

	def derivative
		#: Returns polynomial that is derivative of the present one
		p = []
		(@c.size - 1).times{|i|
			p << @c[i] * (order - i)
		}
		return  Polynomial.new(p)
	end

	def calc(x)
		#: Calculate value for given x
		y = 0.0
		@c.size.times{|i|
			y += @c[i] *  x**(order - i)
		}
		return y
	end

	def print(format = "%20.15f")
		#: Verbose printing of the polynomial
		@c.size.times{|i|
			printf("#{format} * x^%d", @c[i], order - i)
			if i == order
				puts
			else
				puts " +"
			end
		}
		return nil
	end

	def to_s(format = "%.3f")
		s = "Polynomial: "
		@c.size.times{|i|
			s += sprintf("#{format} * x^%d", @c[i], order - i)
			s += " + " unless i == order
		}
		return s
	end

	def order
		# Return ordter of the polynomial
		return @c.size - 1
	end

	def roots # => Array of Complex
		#: Return roots (including complex ones) of the polynomial

		# Companion matrix
		cm = Matrix.zero(order)
		# Fill coefficient (divided by the highest one to get the monic form)
		order.times{|i|
			cm[i,order-1] = -@c[order - i] / @c[0]
		}
		# Sub-diagonal of ones
		(order - 1).times{|i|
			cm[i+1,i] = 1.0
		}
		# Diagonalize companion matrix
		vectors,real,imag = cm.eigensystem

		# Eigenvalues are the roots
		ra = []
		real.m.times{|i|
			ra << Complex(real[i,0], imag[i,0])
		}
		return ra
	end

end
