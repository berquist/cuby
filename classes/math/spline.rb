require "classes/algebra/algebra.rb"

# Cubic spline with optional control of the end points
# By default, a spline ends with points with 2nd derivative = 0
# Optionally, the value of first derivative in the end points can be specified instead

class Spline
	attr_accessor :range

	class Segment
		attr_reader :a, :b, :c, :d

		def initialize(a,b,c,d)
			@a = a
			@b = b
			@c = c
			@d = d
		end

		def calc(x)
			return @a * x**3 + @b * x**2 + @c * x + @d
		end
	end
	
	def initialize(values, deriv_in_a = nil, deriv_in_b = nil)
		n_coeff = (values.size - 1) * 4

		n_eq = n_coeff

		# Matrix
		m = Matrix.zero(n_eq, n_coeff)

		# Right side
		my = Matrix.zero(n_eq, 1)

		# Each segment
		row = 0
		(values.size - 1).times{|i|
			offset = 4 * i
			# Values: first point
			x, y = values[i]
			m[row, offset + 0] = x**3
			m[row, offset + 1] = x**2
			m[row, offset + 2] = x
			m[row, offset + 3] = 1.0
			my[row,0] = y
			row += 1
			# Values: second point
			x, y = values[i + 1]
			m[row, offset + 0] = x**3
			m[row, offset + 1] = x**2
			m[row, offset + 2] = x
			m[row, offset + 3] = 1.0
			my[row,0] = y
			row += 1
		}

		# Midpoints: derivative equality
		(1 .. values.size - 2).each{|i|
			x, y = values[i]

			# First derivative
			offset = (i-1) * 4
			m[row, offset + 0] = 3.0 * x**2
			m[row, offset + 1] = 2.0 * x
			m[row, offset + 2] = 1.0
			offset = (i) * 4
			m[row, offset + 0] = -3.0 * x**2
			m[row, offset + 1] = -2.0 * x
			m[row, offset + 2] = -1.0

			my[row,0] = 0
			row += 1

			# Second derivative
			offset = (i-1) * 4
			m[row, offset + 0] = 6.0 * x
			m[row, offset + 1] = 2.0 
			offset = (i) * 4
			m[row, offset + 0] = -6.0 * x
			m[row, offset + 1] = -2.0

			my[row,0] = 0
			row += 1
		}

		# First point
		x, y = values[0]
		offset = 0
		if deriv_in_a
			# First derivative known
			m[row, offset + 0] = 3.0 * x**2
			m[row, offset + 1] = 2.0 * x
			m[row, offset + 2] = 1.0
			my[row,0] = deriv_in_a
		else
			# Otherwise second set to zero
			m[row, offset + 0] = -6.0 * x
			m[row, offset + 1] = -2.0
			my[row,0] = 0.0
		end
		row += 1

		# Last point
		x, y = values.last
		offset = (values.size - 2) * 4
		if deriv_in_b
			# First derivative known
			m[row, offset + 0] = 3.0 * x**2
			m[row, offset + 1] = 2.0 * x
			m[row, offset + 2] = 1.0
			my[row,0] = deriv_in_b
		else
			# Otherwise second set to zero
			m[row, offset + 0] = -6.0 * x
			m[row, offset + 1] = -2.0
			my[row,0] = 0.0
		end
		row += 1

		# Solve
		poly = (m.inverse * my)

		# Save segments
		@segments = {}

		(values.size - 1).times{|i|
			offset = 4 * i
			# Values: first point
			seg = Segment.new(
				poly[offset + 0, 0],
				poly[offset + 1, 0],
				poly[offset + 2, 0],
				poly[offset + 3, 0]
			)
			@segments[values[i][0] ... values[i+1][0]] = seg
		}

		@range = (values.first[0] .. values.last[0])
	end

	def calculate(x)
		# General range check
		###

		@segments.each{|range, poly|
			if range.cover?(x)
				return poly.calc(x)
			end	
		}
		return nil
	end

end

__END__
# Example:

include Algebra

spline = Spline.new(
	[[0.0,	0.0],
	 [2.95, -1.0],
	 [5.5,	0.0]],
	 nil, 0.0
)


x = 0.0
while x <= 5.5 do
	puts "#{x}\t#{spline.calculate(x)}"
	x += 0.02
end
