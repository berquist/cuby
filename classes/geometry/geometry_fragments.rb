
# Extension of Geometry for manipulation with fragments of the geometry
# e.g. molecular fragments can be obtained from Connectivity.molecules

module GeometryFragments

	def closest_contact(atomlist_a, atomlist_b)
		#: returns distance, and indexes(in geometry) of the two closest atoms
		minimum = 1.0e20
		a = nil
		b = nil
		atomlist_a.each{|i|
			atomlist_b.each{|j|
				raise "Searching for closest contact between fragments, but they overlap" if i == j
				dist = at(i).distance(at(j))
				if dist < minimum
					a = i
					b = j
					minimum = dist
				end
			}
		}
		return [minimum, a, b]
	end

	def closest_contact_geo2(geo2)
		#: returns distance, and indexes(in their geometry) of the two closest atoms
		minimum = 1.0e20
		a = nil
		b = nil
		each_index{|i|
			geo2.each_index{|j|
				dist = at(i).distance(geo2.at(j))
				if dist < minimum
					a = i
					b = j
					minimum = dist
				end
			}
		}
		return [minimum, a, b]
	end

	def find_closest_molecule(molecules, starting_molecule_index)
		#: returns distance, molecule index, indexes of the closest atoms
		raise "more than one molecule must be present" if molecules.size <= 1
		minimum = 1.0e20
		mol = nil
		a = nil
		b = nil
		molecules.each_index{|i|
			next if i == starting_molecule_index
			dist, xa, xb = closest_contact(molecules[starting_molecule_index], molecules[i])
			if dist < minimum
				minimum = dist
				mol = i
				a = xa
				b = xb
			end
		}
		return [minimum, mol, a, b]
	end

end
