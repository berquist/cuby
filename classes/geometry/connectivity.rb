require "classes/algebra/algebra.rb"
require "classes/geometry/geo_grid_cubic.rb"

class Connectivity < SparseMatrix
	# The connectivity is stored as a sparse, symmetric matrix A where
	# bond between atoms i and j is represented as A[i,j] = 1 and
	# A[j,i] = 1

	def initialize(geometry, tolerance = 1.1, force_slow = false)
		# Initialize the matrix, using symmetry
		if geometry.class == Geometry
			super(geometry.size, geometry.size)
		elsif geometry.class == Fixnum
			super(geometry,geometry)
			return
		else
			raise(TypeError, "Connectivity can be build only from geometry")
		end


		# Precalculate distance thresholds
		elements = geometry.elements_in_system
		rcov = {}
		rc_max = 0.0
		elements.each{|e1|
			rcov[e1] = {}
			elements.each{|e2|
				r =(PeriodicTable.covalent_radius(e1) + PeriodicTable.covalent_radius(e2)) * tolerance
				rcov[e1][e2] = r
				rc_max = r if r > rc_max
			}
		}

		if geometry.size < 100 || force_slow
			# Small molecules - simple search
			at_i = nil
			at_j = nil
			geometry.size.times{|i|
				at_i = geometry.at(i)
				i.times{|j|
					at_j = geometry.at(j)
					if at_i.distance(at_j) < rcov[at_i.element][at_j.element]
						# Consider it bond
						self[i,j] = 1.00
						self[j,i] = 1.00
					end
				}
			}
		else
			# Large systems - linear scaling algorithm dividing the molecule
			# using a cubic grid

			# Build grid
			grid = GeoGridCubic.new(geometry, rc_max + 1e-5)

			# Pairs within each cube
			grid.each{|cube|
				cube.size.times{|i|
					at_i = geometry.at(cube[i])
					i.times{|j|
						at_j = geometry.at(cube[j])
						if at_i.distance(at_j) < rcov[at_i.element][at_j.element]
							# Consider it bond
							self[cube[i],cube[j]] = 1.00
							self[cube[j],cube[i]] = 1.00
						end
					}
				}
			}

			# Pairs in adjacent cubes
			grid.each_adjacent_pair{|cube1, cube2|
				cube1.each{|i|
					at_i = geometry.at(i)
					cube2.each{|j|
						at_j = geometry.at(j)
						if at_i.distance(at_j) < rcov[at_i.element][at_j.element]
							# Consider it bond
							self[i,j] = 1.00
							self[j,i] = 1.00
						end
					}
				}
			}
		end
	end

	def Connectivity.resized(conn_old, size_change)
		raise "Connectivity size change must be >= 0" unless size_change >= 0
		newsize = conn_old.m + size_change
		conn_new = Connectivity.new(newsize)
		# Copy the data
		conn_old.each_nonzero_index{|i,j|
			conn_new[i,j] = 1.0
		}

		return conn_new
	end 

	def molecules # => Array of Atomlists
		#: Return list of separate molecules (not connected by covalent bonds)
		#: as array of index lists.

		molecules = []
		frontier_new = [0]
		try = true

		while try
			# go over bonds to find separate molecule
			molecule = []
			while frontier_new.size > 0
				frontier = frontier_new
				frontier_new = []
				frontier.each {|fr|
					molecule << fr unless molecule.include?(fr)
					row_each_nonzero_index(fr){|bound|
						frontier_new << bound unless molecule.include?(bound)
					}
				}
			end

			# save the molecule
			molecules << molecule

			# find another free atom to start at
			try = false
			m.times{|i|
				in_molecule = false
				molecules.each {|mol|
					if mol.include?(i)
						in_molecule = true 
						break
					end
				}
				unless in_molecule
					try = true
					frontier_new = [i]
					break
				end
			}
		end

		# sort atoms in molecules by index
		molecules.each {|mol|
			mol.sort!
		}

		return molecules
	end

	# Some methods are just aliases of the SparseMatrix methods:

	# def hybridization(index) # => Fixnum
	alias :hybridization :row_nonzero_count
	# def bound_atoms(index) # => Atomlist
	alias :bound_atoms :row_nonzero_indices_array
	# def is_bond?(i, j) => Boolean
	alias :is_bond? :nonzero?

	#=======================================================================
	# Writing
	#=======================================================================
	
	def add_bond(i,j)
		self[i,j] = 1.0
		self[j,i] = 1.0
		return nil
	end

	#=======================================================================
	# Iterators
	#=======================================================================

	def each_bond
		each_nonzero_index{|i,j|
			yield(i,j) if j < i
		}
		return nil
	end

	def each_angle
		each_bond{|i,j|
			row_each_nonzero_index(j){|k|
				next if k >= i
				yield(i,j,k)
			}
			row_each_nonzero_index(i){|k|
				next if k >= j
				yield(j,i,k)
			}
		}
		return nil
	end

	def each_torsion
		each_bond{|j,k|
			row_each_nonzero_index(j){|i|
				next if i == k
				row_each_nonzero_index(k){|l|
					next if l == j
					next if i == l
					yield(i,j,k,l)
				}
			}
		}
		return nil
	end

	def one_torsion_per_bond
		each_bond{|j,k|
			found = false
			row_each_nonzero_index(j){|i|
				next if i == k
				row_each_nonzero_index(k){|l|
					next if l == j
					next if i == l
					yield(i,j,k,l)
					found = true
					break if found
				}
				break if found
			}
		}
		return nil
	end

	#!# Inefficient implementation
	# For heavy use, list (SparseMatrix) of pairs should be calculated once and then used
	# Storing excluded pairs would be more efficient in larger systems (1-2) are already stored as bonds,
	# only 1-3 would have to be stored
	def each_nonbonded_above_1_3
		m.times{|i|
			 i.times {|j|
				 next if is_bond?(i,j) # skip bonds
				 next if (row_nonzero_indices_array(i) & row_nonzero_indices_array(j)).size != 0 # no common bonded atoms
				 yield(i,j)
			 }
		}
		return nil
	end

	def each_nonbonded_above_1_4
		each_nonbonded_above_1_3{|i,j|
			is_1_4 = false
			row_nonzero_indices_array(i).each{|k|
				if (row_nonzero_indices_array(k) & row_nonzero_indices_array(j)).size != 0
					is_1_4 = true 
					break
				end
			}
			next if is_1_4
			yield(i,j)
		}
		return nil
	end

	#=======================================================================
	# Lists
	#=======================================================================
	
	def bond_list
		# Array of [i,j] arrays
		bondlist = []
		each_bond{|i,j|
			bondlist << [i,j]
		}
		return bondlist
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def to_s(bond_format="%d-%d", bond_separator=", ")
		#: List of the bonds as a string, atom numbering starts at 1
		s = []
		each_bond{|i,j| s << sprintf(bond_format, i+1, j+1)}
		return s.join(bond_separator)
	end

	def inspect
		return "\nConnectivity:\n" + self.to_s
	end

	#=======================================================================
	# Noncovalent interactions
	#=======================================================================
	
	def each_h_bond(geometry, distance_cutoff = 2.5, angle_cutoff = 60.0, h_bonding_elements = [:O, :N])
		#: Yields triplets of atom indices: donor, hydrogen and acceptor
		m.times{|ai| # traverse acceptors
			next unless h_bonding_elements.include?(geometry.at(ai).element)
			acceptor = geometry.at(ai)
			m.times{|hi| # traverse hydrogens
				next unless geometry.at(hi).element == :H
				hydrogen = geometry.at(hi)
				next if is_bond?(ai,hi) # skip bound pairs
				next if (row_nonzero_indices_array(ai) & row_nonzero_indices_array(hi)).size != 0 # skip 1-3 interactions

				# donor atom
				Cuby::warning("Hydrogen with more than one covalent bonds encountered in H-bond search") if bound_atoms(hi).size > 1
				Cuby::warning("Hydrogen with no covalent bonds encountered in H-bond search") if bound_atoms(hi).size == 0
				di = bound_atoms(hi)[0]
				donor = geometry.at(di)
				next unless h_bonding_elements.include?(donor.element)

				# Cutoffs
				next if acceptor.distance(hydrogen) > distance_cutoff
				angle_d = 180.0 - Coordinate.angle(donor, hydrogen, acceptor) * 180.0 / Math::PI
				next if angle_d > angle_cutoff

				# Passed all tests
				yield(di, hi, ai)
			}
		}
	end

	def h_bond_list(geometry, distance_cutoff = 2.5, angle_cutoff = 60.0, h_bonding_elements = [:O, :N])
		list = []
		each_h_bond(geometry, distance_cutoff, angle_cutoff, h_bonding_elements) {|di, hi, ai|
			list << [di, hi, ai]
		}
		return list
	end

	#=======================================================================
	# Atom types
	#=======================================================================
	
	def atom_types_simple(geometry)
		# Build arrays of elements around
		element_arrays = []
		geometry.each_index{|i|
			# Each atom
			element_arrays[i] = []
			row_each_nonzero_index(i){|j|
				# Each bonded atom
				element_arrays[i] << geometry[j].element
			}
			# Sort the element array
			element_arrays[i].sort!{|a,b|  PeriodicTable.proton_number(a) <=> PeriodicTable.proton_number(b)}
			# Convert the elements to strings
			element_arrays[i].map!{|s| s.to_s }
			# Count occurences
			new = element_arrays[i].uniq
			new.each_index{|j|
				c = 0
				element_arrays[i].each{|s|
					c += 1 if s == new[j]
				}
				new[j] += c.to_s unless c == 1
			}
			element_arrays[i] = new
		}


		atomtypes = []
		geometry.each_index{|i|
			atomtypes[i] = geometry[i].element.to_s + ":" + element_arrays[i].join('')
		}

		return atomtypes
	end

end
