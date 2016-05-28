class Polarizability
	attr_reader :matrix

	def initialize(units)
		@matrix = Matrix.zero(3)
		@units = units
	end

	def to_s
		s = "Polarizability (#{@units}):\n"
		s += @matrix.to_s
		s += "\n   1/3*trace: #{'%.5f' % (@matrix.trace / 3)}"
		return s
	end
end
