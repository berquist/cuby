require "classes/objects/molecular_orbital.rb"

class MolecularOrbitals < Array

	def homo
		# Search for HOMO and return it
		max = first.energy_au
		orb = first
		each{|mo|
			if mo.energy_au > max && mo.occupation > 0
				max = mo.energy_au
				orb = mo
			end
		}
		return orb
	end

	def lumo
		# Search for LUMO and return it
		min = last.energy_au
		orb = last
		each{|mo|
			if mo.energy_au < min && mo.occupation == 0
				min = mo.energy_au
				orb = mo
			end
		}
		return orb
	end

	def sort_by_e!
		sort!{|a,b| a.energy_au <=> b.energy_au}
	end
end
