################################################################################
#
# Class InternalCoordinate
#
# Author: Jan Rezac
# Date created: 2008-11-27
# License: Cuby license
# Description: Classes for work with internal coordinates
# Status: Works
#
################################################################################

#: Three classes for representation of internal coordinate:
#* InternalCoordinate: works directly on Coordinates or Atoms
#* GroupsInternalCoord: Works on groups of Coordinates or Atoms, groups are represented by their center of mass for the calculations
#* InternalCoordinateI: works with atom indexes in specified Geometry

class InternalCoordError < RuntimeError
end

class InternalCoordinate

	### Add options:
	@@check_identical_atoms = true
	@@report_undefined_torsions = false

	attr_accessor 	:atoms
	attr_reader	:type
	attr_accessor :name # The coordinate can be named, but this feature is not used by default

	def initialize(type, *atoms)
		#: Types:
		#* :distance
		#* :distance_inv
		#* :angle
		#* :torsion
		#* :distance_difference
		@type = type
		case @type
		when :distance
			raise(InternalCoordError, "Two cartesians are required for calculation of distance") unless atoms.size == 2
		when :distance_inv
			raise(InternalCoordError, "Two cartesians are required for calculation of inverse distance") unless atoms.size == 2
		when :angle
			raise(InternalCoordError, "Three cartesians are required for calculation of angle") unless atoms.size == 3
		when :torsion
			raise(InternalCoordError, "Four cartesians are required for calculation of torsion angle") unless atoms.size == 4
		when :distance_difference
			raise(InternalCoordError, "Four cartesians are required for calculation of distance difference") unless atoms.size == 4
		else
			raise(InternalCoordError, "Unknown type of internal coordinate")
		end

		if @@check_identical_atoms && @type != :distance_difference
			atoms.each_index{|i|
				i.times{|j|
					if atoms[i].object_id == atoms[j].object_id
						raise(InternalCoordError, "Wrong definition of internal coordinate (#{@type}) - identical atoms #{i} and #{j}")
					end
				}
			}
		end

		@atoms = atoms
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def size # => integer
		#: Get dimension of the coordinate
		return @atoms.size
	end

	def InternalCoordinate.size_of(type)
		#: Get dimension of a coordinate of specified type
		case type
		when :distance
			return 2
		when :distance_inv
			return 2
		when :angle
			return 3
		when :torsion
			return 4
		when :distance_difference
			return 4
		else
			raise(InternalCoordError, "Unknown type of internal coordinate")
		end
	end

	#=======================================================================
	# Calculations
	#=======================================================================

	def value # => float
		#: Returns value of the internal coordinate
		case @type
		when :distance
			return @atoms[0].distance(@atoms[1])
		when :distance_inv
			return 1.0 / @atoms[0].distance(@atoms[1])
		when :angle
			return Coordinate.angle(*@atoms)
		when :torsion
			return Coordinate.torsion(*@atoms)
		when :distance_difference
			return @atoms[0].distance(@atoms[1]) - @atoms[2].distance(@atoms[3])
		end
	end

	def grad_on_atoms(grad) # => Array of Coordinates
		#: Returns cartesian gradients on the atoms calculated from gradient in the internal coordinate
		gradients = @atoms.map{|list| Coordinate.new }
		case @type
		when :distance
			i,j = @atoms
			u = i - j
			gradients[0] = grad / u.abs * u
			gradients[1] = -gradients[0]
		when :angle
			i, j, k = @atoms

			u = i - j
			v = k - j

			ru = u.abs
			rv = v.abs
			dot = u.dot(v)

			gradients[0] = grad / (1.0 - dot**2 / ru**2 / rv**2)**0.5 * -(v/ru/rv - u*dot/ru**3/rv)
			gradients[2] = grad / (1.0 - dot**2 / ru**2 / rv**2)**0.5 * -(u/ru/rv - v*dot/ru/rv**3)
			gradients[1] = -gradients[0] - gradients[2]
		when :torsion
			i, j, k, l = @atoms
			u = j - i
			v = k - j
			w = l - k
			t = i - k
			s = j - l
			v_abs = v.abs

			# zero distance check
			[u,v,w,t,s].each{|vec|
				return gradients if vec.abs == 0.0
			}

			# linear case check
			a1 = u.angle(v)
			a2 = v.angle(w)
			if a1 == 0.0 || a1 == Math::PI ||
			   a2 == 0.0 || a2 == Math::PI
				return gradients
			end

			cuv = u.cross_product(v)
			cvw = v.cross_product(w)

			ex1 = v_abs * u.dot(cvw) 
			ex2 = cuv.dot(cvw) 

			denom = (ex1**2+ex2**2)

			gradients[0] = 
				v.cross_product(cvw)*ex1 - 
				cvw*v_abs*ex2
			gradients[3] = 
				v.cross_product(cuv)*ex1 + 
				cuv*v_abs*ex2

			gradients[1] = -u.dot(cvw)/v_abs*ex2*v + 
				cvw*v_abs*ex2 + 
				u.cross_product(w)*v_abs*ex2 +
				t.cross_product(cvw)*ex1 - 
				cuv.cross_product(w)*ex1

			gradients[2] = u.dot(cvw)/v_abs*ex2*v + 
				u.cross_product(s)*v_abs*ex2 +
				u.cross_product(cvw)*ex1 - 
				cuv.cross_product(s)*ex1

			4.times {|i| gradients[i] *= grad / denom}
		when :distance_difference
			# distance 1
			u = @atoms[0] - @atoms[1]
			gradients[0] = grad / u.abs * u
			gradients[1] = -gradients[0]
			v = @atoms[2] - @atoms[3]
			gradients[2] = -grad / v.abs * v
			gradients[3] = -gradients[2]
		end

		return gradients
	end

	def grad_to_atoms(grad)
		#: gradient in the internal coordinate is applied to atoms in the geometry
		g = grad_on_atoms(grad)
		@atoms.each_index{|i|
			@atoms[i].grad += g[i]
		}
		return nil
	end

	def grad_to_grad_array(grad, geometry, array)
		#: gradient in the internal coordinate is applied to array of gradients
		g = grad_on_atoms(grad)
		@atoms.each_index{|i|
			index = geometry.index(@atoms[i])
			raise(InternalCoordError, "Atom in the internal coordinate is not a member of the geometry") if index == nil
			array[index] += g[i]
		}
		return nil
	end
end

class GroupsInternalCoord

	attr_reader :lists
	attr_accessor :name # The coordinate can be named, but this feature is not used by default

	#=======================================================================
	# Initialize
	#=======================================================================

	def initialize(type, geometry, *lists)
		#: Types:
		#* :distance
		#* :distance_inv
		#* :angle
		#* :torsion
		@type = type
		case @type
		when :distance
			raise(InternalCoordError, "Two cartesians are required for calculation of distance") unless lists.size == 2
		when :distance_inv
			raise(InternalCoordError, "Two cartesians are required for calculation of inverse distance") unless lists.size == 2
		when :angle
			raise(InternalCoordError, "Three cartesians are required for calculation of angle") unless lists.size == 3
		when :torsion
			raise(InternalCoordError, "Four cartesians are required for calculation of torsion angle") unless lists.size == 4
		when :distance_difference
			raise(InternalCoordError, "Four cartesians are required for calculation of distance difference") unless lists.size == 4
		else
			raise(InternalCoordError, "Unknown type of internal coordinate")
		end
		lists.each {|list|
			raise(InternalCoordError, "List of atoms should be array") unless list.class == Array
		}
		@lists = lists
		@geometry = geometry
		centers = @lists.map{|list| @geometry.geometry_from_list(list).center_of_mass}
		@ic = InternalCoordinate.new(@type, *centers)
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def size # => integer
		#: Get dimension of the coordinate
		return @lists.size
	end

	def type # => Symbol
		#: Get type of the coordinate
		return @ic.type
	end

	#=======================================================================
	# Calculations
	#=======================================================================

	def value # => float
		#: Returns value of the internal coordinate
		centers = @lists.map{|list| @geometry.geometry_from_list(list).center_of_mass}
		@ic.atoms = centers
		return @ic.value
	end

	def grad_on_centers(grad) # => Array of Coordinates
		#: Returns cartesian gradients on the centers calculated from gradient in the internal coordinate
		centers = @lists.map{|list| @geometry.geometry_from_list(list).center_of_mass}
		@ic.atoms = centers
		return @ic.grad_on_atoms(grad)
	end

	def grad_to_geometry(grad)
		#: gradient in the internal coordinate is applied to atoms in the geometry
		g = grad_on_centers(grad)
		@geometry.each_in_list(@lists[i]) {|atom|
			atom.grad += g[i] / @lists[i].size
		}
		return nil
	end

	def grad_to_grad_array(grad, array)
		#: gradient in the internal coordinate is applied to array of gradients
		g = grad_on_centers(grad)
		@lists.each_index {|i|
			@lists[i].each {|index|
				array[index] += g[i] / @lists[i].size
			}
		}
		return nil
	end

end

class InternalCoordinateI
	#: Here the atoms are specified by their index in geometry

	attr_reader :indexes
	attr_accessor :name # The coordinate can be named, but this feature is not used by default

	#=======================================================================
	# Initialize
	#=======================================================================

	def initialize(type, geometry, *indexes)
		#: Types:
		#* :distance
		#* :distance_inv
		#* :angle
		#* :torsion
		@type = type
		case @type
		when :distance
			raise(InternalCoordError, "Two cartesians are required for calculation of distance") unless indexes.size == 2
		when :distance_inv
			raise(InternalCoordError, "Two cartesians are required for calculation of inverse distance") unless indexes.size == 2
		when :angle
			raise(InternalCoordError, "Three cartesians are required for calculation of angle") unless indexes.size == 3
		when :torsion
			raise(InternalCoordError, "Four cartesians are required for calculation of torsion angle") unless indexes.size == 4
		when :distance_difference
			raise(InternalCoordError, "Four cartesians are required for calculation of distance difference") unless indexes.size == 4
		else
			raise(InternalCoordError, "Unknown type of internal coordinate")
		end
		@indexes = indexes
		@geometry = geometry
		atoms = @indexes.map{|i| @geometry[i]}
		@ic = InternalCoordinate.new(@type, *atoms)
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def size # => integer
		#: Get dimension of the coordinate
		return @indexes.size
	end

	def type # => Symbol
		#: Get type of the coordinate
		return @ic.type
	end

	def atoms # => Array of Atom
		#: Get array of the atoms
		atoms = @indexes.map{|i| @geometry[i]}
		return atoms
	end

	#=======================================================================
	# Calculations
	#=======================================================================

	def value # => float
		#: Returns value of the internal coordinate
		return @ic.value
	end

	def grad_on_atoms(grad) # => Array of Coordinates
		#: Returns cartesian gradients on the centers calculated from gradient in the internal coordinate
		return @ic.grad_on_atoms(grad)
	end

end
