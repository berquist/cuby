# Classes used for the sorage of the isotope data in the yaml file

class Isotope
	attr_reader :number, :mass, :abundance

	def mass=(mass)
		@mass = mass
	end
end

class IsotopeList < Hash
end

