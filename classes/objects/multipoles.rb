class Multipoles < Hash

	def initialize
	end

	class Dipole
		attr_accessor :value # Coordinate type
		attr_accessor :derivative # Array of XYZ coordinates

		def initialize(x,y,z)
			@value = Coordinate[x,y,z]
		end

		def to_s(unit = :debye)
			case unit
			when :debye
				v = @value * CUBY2DEBYE
				unit = "Debye"
			else
				v = @value
				unit = "A*e"
			end
			s = sprintf("Dipole: [%.5f, %.5f, %.5f] %s", v.x, v.y, v.z, unit)
			s += "\n"
			s += sprintf("Dipole size: %.5f %s", v.abs, unit)
			return s
		end

		def size
			return @value.abs
		end

		def [](index)
			return @value[index]
		end
	end
end
