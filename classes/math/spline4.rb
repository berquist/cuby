require "classes/algebra/algebra.rb"

# Fourth-order spline
# Compared to cubic spline, one additional condition can be applied to the end points
# Three out of four possible conditions must be set:
# :start_d1 - first derivative in the first point
# :start_d2 - second derivative in the first point
# :end_d1 - first derivative in the last point
# :end_d2 - second derivative in the last point

class Spline4
	attr_accessor :range

	COND_KEYS = [:start_d1, :start_d2, :end_d1, :end_d2]

	class Segment
		attr_reader :a, :b, :c, :d, :e

		def initialize(a,b,c,d,e)
			@a = a
			@b = b
			@c = c
			@d = d
			@e = e
		end

		def calc(x)
			return @a * x**4 + @b * x**3 + @c * x**2 + @d * x + @e
		end
	end
	
	POLY_S = 5

	def initialize(values, conditions)
		# Check the conditions
		raise "Three conditions must be specified" unless conditions.size == 3
		raise "Invalid condition name(s): #{(conditions.keys - COND_KEYS).join{", "}}" if (conditions.keys - COND_KEYS).size > 0

		n_coeff = (values.size - 1) * POLY_S

		n_eq = n_coeff

		# Matrix
		m = Matrix.zero(n_eq, n_coeff)

		# Right side
		my = Matrix.zero(n_eq, 1)

		# Each segment
		row = 0
		(values.size - 1).times{|i|
			offset = POLY_S * i
			# Values: first point
			x, y = values[i]
			m[row, offset + 0] = x**4
			m[row, offset + 1] = x**3
			m[row, offset + 2] = x**2
			m[row, offset + 3] = x
			m[row, offset + 4] = 1.0
			my[row,0] = y
			row += 1
			# Values: second point
			x, y = values[i + 1]
			m[row, offset + 0] = x**4
			m[row, offset + 1] = x**3
			m[row, offset + 2] = x**2
			m[row, offset + 3] = x
			m[row, offset + 4] = 1.0
			my[row,0] = y
			row += 1
		}

		# Midpoints: derivative equality
		(1 .. values.size - 2).each{|i|
			x, y = values[i]

			# First derivative
			offset = (i-1) * POLY_S
			m[row, offset + 0] = 4.0 * x**3
			m[row, offset + 1] = 3.0 * x**2
			m[row, offset + 2] = 2.0 * x
			m[row, offset + 3] = 1.0
			offset = (i) * POLY_S
			m[row, offset + 0] = -4.0 * x**3
			m[row, offset + 1] = -3.0 * x**2
			m[row, offset + 2] = -2.0 * x
			m[row, offset + 3] = -1.0

			my[row,0] = 0
			row += 1

			# Second derivative
			offset = (i-1) * POLY_S
			m[row, offset + 0] = 12.0 * x**2
			m[row, offset + 1] = 6.0 * x
			m[row, offset + 2] = 2.0 
			offset = (i) * POLY_S
			m[row, offset + 0] = -12.0 * x**2
			m[row, offset + 1] = -6.0 * x
			m[row, offset + 2] = -2.0

			my[row,0] = 0
			row += 1

			# Third derivative
			offset = (i-1) * POLY_S
			m[row, offset + 0] = 24.0 * x
			m[row, offset + 1] = 6.0 
			offset = (i) * POLY_S
			m[row, offset + 0] = -24.0 * x
			m[row, offset + 1] = -6.0

			my[row,0] = 0
			row += 1
		}

		# First point
		x, y = values[0]
		offset = 0
		if conditions[:start_d1]
			# First derivative known
			m[row, offset + 0] = 4.0 * x**3
			m[row, offset + 1] = 3.0 * x**2
			m[row, offset + 2] = 2.0 * x
			m[row, offset + 3] = 1.0
			my[row,0] = conditions[:start_d1]
			row += 1
		end
		if conditions[:start_d2]
			# Second derivative
			m[row, offset + 0] = -12.0 * x**2
			m[row, offset + 1] = -6.0 * x
			m[row, offset + 2] = -2.0 
			my[row,0] = conditions[:start_d2]
			row += 1
		end

		# Last point
		x, y = values.last
		offset = (values.size - 2) * POLY_S
		if conditions[:end_d1]
			# First derivative known
			m[row, offset + 0] = 4.0 * x**3
			m[row, offset + 1] = 3.0 * x**2
			m[row, offset + 2] = 2.0 * x
			m[row, offset + 3] = 1.0
			my[row,0] = conditions[:end_d1]
			row += 1
		end
		if conditions[:end_d2]
			# Second derivative
			m[row, offset + 0] = -12.0 * x**2
			m[row, offset + 1] = -6.0 * x
			m[row, offset + 2] = -2.0 
			my[row,0] = conditions[:end_d2]
			row += 1
		end

		# Solve
		poly = (m.inverse * my)

		# Save segments
		@segments = {}

		(values.size - 1).times{|i|
			offset = POLY_S * i
			# Values: first point
			seg = Segment.new(
				poly[offset + 0, 0],
				poly[offset + 1, 0],
				poly[offset + 2, 0],
				poly[offset + 3, 0],
				poly[offset + 4, 0]
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

	def calc(x)
		calculate(x)
	end
end

__END__
Example:

include Algebra

spline = Spline4.new(
	[[0.0,	0.0],
	 [2.95, -1.0],
	 [5.5,	0.0]],
	 {:start_d1 => 0, :end_d1 => 0, :end_d2 => 0}
)


x = 0.0
while x <= 5.5 do
	puts "#{x}\t#{spline.calculate(x)}"
	x += 0.02
end
