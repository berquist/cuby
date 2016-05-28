require "classes/algebra/algebra.rb"

# Cubic Hermite spline, with optional monotonicity fix
# Algorithm from https://en.wikipedia.org/wiki/Monotone_cubic_interpolation

class SplineHermite
	attr_accessor :range

	class Segment
		attr_reader :x0, :x1, :y0, :y1, :m0, :m1

		def initialize(x0, x1, y0, y1, m0, m1)
			@x0, @x1, @y0, @y1, @m0, @m1 = x0, x1, y0, y1, m0, m1
		end

		def calc(x)
			t = (x - @x0) / (@x1 - @x0)
			# Hermitre base functions
			h00 = (1.0 + 2.0*t)*(1.0 - t)**2
			h10 = t * (1.0 - t)**2
			h01 = t**2 * (3.0 - 2.0*t)
			h11 = t**2 * (t - 1.0)

			return h00 * @y0 + h10*(@x1-@x0)*@m0 + h01*@y1 + h11*(@x1-@x0)*@m1
		end
	end
	
	def initialize(values, monotone = false, start_deriv = nil, end_deriv = nil)
		@range = (values.first[0] .. values.last[0])
		n_seg = values.size - 1

		# X and Y arrays
		x = []
		y = []
		values.each{|vx,vy|
			x << vx
			y << vy
		}

		# Secant in all points but end
		delta = []
		for i in (0..(values.size - 2)) do
			delta[i] = (y[i+1] - y[i]) / (x[i+1]-x[i])
		end

		# Tangents as average of the secants
		m = []
		for i in (1..(values.size - 2)) do
			m[i] = (delta[i-1] + delta[i]) / 2
		end
		# End points - take linear slope
		m[0] = delta[0]
		m[values.size - 1] = delta[values.size - 2]
		# Or set the derivative manually
		m[0] = start_deriv if start_deriv
		m[values.size - 1] = end_deriv if end_deriv
		
		# Ensure monotonous
		if monotone
			for i in (0..(values.size - 2)) do
				if y[i] == y[i+1]
					m[i] = m[i+1] = 0.0
				else
					a = m[i]/delta[i]
					b = m[i+1]/delta[i]
					m[i] = 3.0 * delta[i] if a > 3
					m[i+1] = 3.0 * delta[i] if b > 3
				end
			end
		end

		# Save segments
		@segments = {}
		n_seg.times{|i|
			seg = Segment.new(values[i][0], values[i+1][0], values[i][1], values[i+1][1], m[i], m[i+1])
			@segments[values[i][0] ... values[i+1][0]] = seg
		}

	end

	def calculate(x)
		# General range check
		unless @range.cover?(x)
			raise "x outsides of range on which the segmented line is defined"
		end

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
require "classes/math/spline"

include Algebra

valueset = [
	[0.0,	0.0],
	[1.0, 0.5],
	[2.0,	2.0],
	[3.0,	0.5],
#	[4.0,	0.0],
	[5.5,	0.0]
]

spline = Spline.new(valueset)
spline2 = SplineHermite.new(valueset, false, 0.0, 0.0)
spline3 = SplineHermite.new(valueset, true, 0.0, 0.0)


x = 0.0
while x <= 5.5 do
	puts "#{x}\t#{spline.calculate(x)}\t#{spline2.calculate(x)}\t#{spline3.calculate(x)}"
	x += 0.02
end
