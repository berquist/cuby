class MolecularOrbital
	attr_reader :name
	attr_reader :energy_au
	attr_reader :occupation

	def initialize(name, energy_au, occupation)
		@name = name
		@energy_au = energy_au
		@occupation = occupation
	end

	def to_s
		sprintf("%10s%16.6f a.u. %10.3f", @name, @energy_au, @occupation)
	end
end
