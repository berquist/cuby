
# Segmented line - linear interpolation between predefined points
# 1:1 replacement for splines

class SegmentedLine
	attr_accessor :range

	class Segment
		attr_reader :a, :b

		def initialize(a, b)
			@a = a
			@b = b
		end

		def calc(x)
			return @a * x + @b
		end
	end
	
	def initialize(values)
		@range = (values.first[0] .. values.last[0])
		n_seg = values.size - 1

		@segments = {}

		n_seg.times{|i|
			a = (values[i+1][1] - values[i][1]) / (values[i+1][0] - values[i][0])
			b = values[i][1]
			seg = Segment.new(a,b)
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
				return poly.calc(x - range.begin)
			end	
		}
		return nil
	end

end

__END__
# Example:

arr = []
x = 0.0
while x < 5.5 do
	arr << [x, Math::sin(x)]
	x += 0.5
end
arr << [5.5, 0.0]

spline = SegmentedLine.new(arr)


x = 0.0
while x <= 5.5 do
	puts "#{x}\t#{spline.calculate(x)}"
	x += 0.02
end
