
class ContinuousValence
	attr_reader :valence
	attr_reader :gradient
	attr_reader :atoms

	def initialize(geometry, do_gradient = false, atoms_in_range = false)
		#: Continuous valence is a smooth function that indicates how many atoms
		#: are within covalent distance. List of the contributing atoms and derivatives
		#: of their contribution to the resulting number are stored if requested.

		# Initialize data
		@valence = Array.new(geometry.size, 0.0)
		if do_gradient
			@gradient = []
			@valence.size.times{|i|
				@gradient[i] = {}
			}
		end
		if atoms_in_range
			@atoms = []
			@valence.size.times{|i|
				@atoms[i] = []
			}
		end

		geometry.each_index{|i|
			atom_i = geometry.at(i)
			r_i = PeriodicTable.covalent_radius(atom_i.element)
			i.times{|j|
				atom_j = geometry.at(j)
				r_j = PeriodicTable.covalent_radius(atom_j.element)
				r0 = r_i + r_j
				r1 = r0 * 1.6

				r = atom_i.distance(atom_j)
				if r > r1
					next # No contribution
				end

				# Save contributing atoms
				if atoms_in_range
					@atoms[i] << j
					@atoms[j] << i
				end

				if r < r0 # Covalent distances
					@valence[i] += 1.0
					@valence[j] += 1.0
					next # No gradient
				end

				x = (r - r0) / (r1 - r0)
				fsw = 1.0-(-20.0*x**7 + 70*x**6 -84.0*x**5 + 35.0*x**4)

				# Valence
				@valence[i] += fsw
				@valence[j] += fsw

				# Gradient along distance
				if do_gradient
					fsw_d = -(-140.0*x**6 + 420.0*x**5 - 420.0*x**4 + 140.0*x**3) / (r1 - r0)
					@gradient[i][j] = fsw_d
					@gradient[j][i] = fsw_d
				end
			}
		}
	end

	def contribution(atom_i, atom_j)
		r_i = PeriodicTable.covalent_radius(atom_i.element)
		r_j = PeriodicTable.covalent_radius(atom_j.element)
		r0 = r_i + r_j
		r1 = r0 * 1.6

		r = atom_i.distance(atom_j)
		if r > r1
			return 0.0
		end

		if r < r0 # Covalent distances
			return 1.0
		end

		x = (r - r0) / (r1 - r0)
		return 1.0-(-20.0*x**7 + 70*x**6 -84.0*x**5 + 35.0*x**4)
	end

	def contribution_d(atom_i, atom_j)
		r_i = PeriodicTable.covalent_radius(atom_i.element)
		r_j = PeriodicTable.covalent_radius(atom_j.element)
		r0 = r_i + r_j
		r1 = r0 * 1.6

		r = atom_i.distance(atom_j)
		if r > r1
			return 0.0
		end

		if r < r0 # Covalent distances
			return 0.0
		end

		x = (r - r0) / (r1 - r0)
		return -(-140.0*x**6 + 420.0*x**5 - 420.0*x**4 + 140.0*x**3) / (r1 - r0)
	end


end
