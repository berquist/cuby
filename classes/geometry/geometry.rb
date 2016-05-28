require "yaml"
require 'digest/md5'

require "classes/geometry/atom.rb"

require "classes/geometry/geometry_selections.rb"
require "classes/geometry/geometry_read_write.rb"
require "classes/geometry/geometry_database.rb"
require "classes/geometry/geometry_smiles.rb"
require "classes/geometry/geometry_fragments.rb"
require "classes/geometry/geometry_pdb.rb"
require "classes/geometry/geometry_rmsd.rb"

require "classes/geometry/connectivity.rb"
require "classes/geometry/rings.rb"

class Geometry < Array
	include GeometrySelections
	include GeometryReadWrite
	include GeometryFragments
	include GeometryPdb
	include GeometryRmsd

	attr_accessor :info

	#=======================================================================
	# Constructors
	#=======================================================================
	
	def initialize
		@info = {}
	end
	
	def Geometry.load_from_settings(settings, keyword = :geometry)
		#: If the argument is a valid file name, it reads the geometry from the file,
		#: otherwise it attempts to read geometry from a string.

		string = settings[keyword]

		if string.class == String
			if FileTest.file?(string)
				# Try file first
				geo = Geometry.from_file(string)
			elsif FileTest.file?(File.expand_path(string))
				# File, with expanded path
				geo = Geometry.from_file(File.expand_path(string))
			elsif string =~ /^parent_block/i
				# load from parent block
				Cuby::error("Requested reading geometry from parent block but there is no parent") unless settings.parent
				geo = Geometry.load_from_settings(settings.parent)
			elsif string =~ /^smiles:/i
				# Geometry to be built from SMILES
				smiles = string.gsub(/^smiles:/i,"")
				geo = GeometrySmiles.build(smiles, settings)
			elsif string =~ /[^:]+:[^:]/
				# Geometry from a database
				geo = GeometryDatabase.get(string, settings)
			else
				# String representation of the geometry
				begin
					geo = Geometry.from_s(string)
				rescue
					if File.basename(string) =~ /^[0-9,a-z,A-Z,-,_]+\.[0-9,a-z,A-Z]+$/
						# The string looks like a filename, but the file was not found
						Cuby::error("Error reading geometry, the entry in the input seems to be a filename\nbut the file was not found")
					else
						raise
					end
				end
			end
		else
			raise(TypeError, "Geometry.load argument must be a string")
		end

		if geo.size == 0
			Cuby::error("Error reading geometry: no atoms found")
		end

		if settings.set?(:geometry_update_coordinates)
			# Update the coordinates in the geometry from another geometry
			# Read the second geometry
			geo_crd = Geometry.from_file(File.expand_path(settings[:geometry_update_coordinates]))
			# Check size
			Cuby::error "Size of geometries in 'geometry' and 'geometry_update_coordinates' must be the same." unless geo.size == geo_crd.size
			# Check matching atoms
			ok = true
			geo.each_index{|i| ok &= geo[i].element == geo_crd[i].element}
			Cuby::recommendation "Elements in 'geometry_update_coordinates' do not match those in 'geometry'" unless ok
			# Update
			geo.copy_coordinates_from(geo_crd)
		end

		if settings.set?(:geometry_rotate)
			roatations = settings[:geometry_rotate].downcase.strip.split(/\s*,\s*/)
			Cuby::error "geometry_rotate expects three comma separated values" unless roatations.size == 3
			roatations.each_index{|i|
				angle = (roatations[i] == 'rand') ? rand : roatations[i].to_f
				axis = Coordinate.new
				axis[i] = 1.0
				geo.rotate_angle_axis!(angle / 180.0 * Math::PI, axis)
			}

		end

		case settings[:geometry_reorder]
		when :reverse
			temp = geo.dup
			geo.clear
			temp.each{|atom| geo.unshift(atom)}
		when :random
			geo.shuffle!
		end

		return geo
	end

	def Geometry.from_file(filename)
		g = Geometry.new
		g.read_file(File.expand_path(filename))
		return g
	end
	

	#=======================================================================
	# Cuby geometry printing/reading
	#=======================================================================

	def to_s
		return map{|atom| atom.to_s}.join("\n")
	end

	def Geometry.from_s(string)
		geo = Geometry.new
		string.each_line{|line|
			next if line =~ /^\s*$/ # ignore blank lines
			next if line =~ /^\s*#/ # ignore comments
			geo << Atom.from_s(line)
		}
		geo.each_index{|i| geo[i].properties[:file_index] = i}
		return geo
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def inspect
		return to_s
	end

	def deep_copy
		geo = Geometry.new
		each{|atom|
			geo << atom.deep_copy
		}
		geo.info = self.info.dup
		return geo
	end

	def checksum
		incr_digest = Digest::MD5.new()
		each{|atom| incr_digest << atom.to_s }
		return incr_digest.hexdigest
	end


	#=======================================================================
	# Geometry as Matrix
	#=======================================================================
	
	def to_matrix
		m = Matrix.zero(self.size, 3)
		each_index{|i|
			atom = at(i)
			m[i,0], m[i,1], m[i,2] = atom.x, atom.y, atom.z
		}
		return m
	end

	#=======================================================================
	# Element and atom info
	#=======================================================================
	
	def elements_in_system
		#: Returns list of elements found in the geometry
		list = []
		each {|atom|
			list |= [atom.element]
		}
		return list
	end

	def count_elements
		#: returns a has of element => count
		list = {}
		each {|atom|
			if list[atom.element]
				list[atom.element] += 1
			else
				list[atom.element] = 1
			end
		}
		return list
	end

	def ghost_atoms?
		#: Returns true when there is a ghost atom in the geometry
		each{|atom|
			return true if atom.properties[:ghost]
		}
		return false
	end

	def dummy_atoms?
		#: Returns true when there is a dummy atom in the geometry
		each{|atom|
			return true if atom.properties[:dummy]
		}
		return false
	end

	def atom_numbering_by_element
		atom_numbers = []
			
		elements_c = {}
		each_index{|i|
			atom = at(i)
			e = atom.element
			if elements_c[e]
				elements_c[e] = elements_c[e] + 1
			else
				elements_c[e] = 1
			end
			atom_numbers[i] =  elements_c[e]
		}
		return atom_numbers
	end


	#=======================================================================
	# YAML support
	#=======================================================================
	
	yaml_as "tag:cuby.molecular.cz,2009:geometry"

	def to_yaml(opts = {})
		YAML::quick_emit( self.object_id, opts ) do |out|
			out.scalar(taguri, self.to_s, :literal )
		end
	end

	YAML::add_domain_type( "cuby.molecular.cz,2009", "geometry") do  |type, val|
		Geometry.from_s(val)
	end

	#=======================================================================
	# Movement
	#=======================================================================
	
	def translate!(coordinate) # => nil
		#: Adds a vector (instance of Coordinate) to coordinates of each atom
		#: in the geometry
		each {|atom|
			atom.plus!(coordinate)
		}
		return nil
	end

	def translate(coordinate) # => Geometry
		#: Returns translated copy of the geomery
		tra = self.deep_copy
		tra.translate!(coordinate)
		return tra
	end

	def rotate_angle_axis!(angle, axis)
		each {|atom|
			atom.rotate_angle_axis!(angle,axis)
		}
		return nil
	end

	def rotate!(matrix)
		each {|atom|
			atom.rotate_by_matrix!(matrix)
		}
		return nil
	end

	def orient!(l1,l2,l3)
		# Align selected three atoms / centers to origin, x axis and x-y plane
		g1 = self.geometry_from_list(l1)
		g2 = self.geometry_from_list(l2)
		g3 = self.geometry_from_list(l3)

		# Translate to origin
		c1 = g1.center
		self.translate!(-c1)

		# Rotate to x axis
		c2 = g2.center
		x_axis = Coordinate[-1,0,0]
		angle = x_axis.angle(c2)
		if angle == Math::PI
			# In axis, but reverse
			# Normal can not be used
			rot_axis = Coordinate[0,1,0]
			rotator = CoordinateVector.rotator_from_angle_axis(-angle,rot_axis)
			self.rotate!(rotator)
		elsif angle > 0.0
			rot_axis = x_axis.cross_product(c2)
			rotator = CoordinateVector.rotator_from_angle_axis(-angle,rot_axis)
			self.rotate!(rotator)
		end

		# Rotate to plane
		if size > 2 # Skip for diatomics
			c1 = g1.center
			c2 = g2.center
			c3 = g3.center
			angle = Coordinate.torsion(Coordinate[0,1,0],c1,c2,c3)
			rotator = CoordinateVector.rotator_from_angle_axis(angle,Coordinate[1,0,0])
			self.rotate!(rotator)
		end

		return nil
	end

	#=======================================================================
	# Geometry operations
	#=======================================================================
	
	def copy_coordinates_from(geometry)
		each_index{|i|
			at(i).set_coord(geometry[i])
		}
	end

	#=======================================================================
	# Centroids
	#=======================================================================

	def weighted_center(weight_type, weights = nil) # => coordinate
		#: Calculate center of the system using weights specified by weight_type parameter:
		#* :unit - equal weights of 1.0
		#* :mass - Mass of atom
		#* :array - weights taken from second parameter, Array of floats

		case weight_type
		when :mass
		when :unit
		when :array
			raise "# of weights != # of atoms" unless weights.size == natom
		else
			raise "Unknown weight type"
		end

		center = Coordinate.new
		totalweight = 0.0
		each_index { |i|
		       atom = at(i)	
			case weight_type
			when :mass
				weight = atom.mass
			when :unit
				weight = 1.0
			when :array
				weight = weights[i]
			end

			center.x += atom.x * weight
			center.y += atom.y * weight
			center.z += atom.z * weight
			totalweight += weight

		}
		center.x /= totalweight
		center.y /= totalweight
		center.z /= totalweight

		return center
	end

	def com
		return weighted_center(:mass)
	end

	def center_of_mass
		return weighted_center(:mass)
	end

	def center
		return weighted_center(:unit)
	end

	def weighted_center_atomlist(atomlist,weight_type)
		# Weighted center of selected atoms
		g = self.geometry_from_list(atomlist)
		return g.weighted_center(weight_type)
	end

	def center_of_mass_atomlist(atomlist)
		return weighted_center_atomlist(atomlist, :mass)
	end

	#=======================================================================
	# Operators
	#=======================================================================
	
	def +(geometry2)
		#: Concatenation of geometry objects
		g = Geometry.new
		each {|atom|
			g << atom
		}
		geometry2.each {|atom|
			g << atom
		}
		return g
	end
	
	#=======================================================================
	# Iterators
	#=======================================================================
	
	def each_with_index
		each_index{|i|
			yield(at(i), i)
		}
	end

	#=======================================================================
	# Vector access
	#=======================================================================
	
	def to_vector # => Vector
		#: Returns Vector containing atomic coordinates in format
		#: (atom0.x, atom0.y, atom0.z, atom1.x, ...)
		array = []
		each {|atom|
			array << atom.x
			array << atom.y
			array << atom.z
		}
		return Vector.from_array(array)
	end

	def update_from_vector!(vector) # => nil
		#: Read atomic coordinates from Vector of size natom * 3 in format
		#: (atom0.x, atom0.y, atom0.z, atom1.x, ...)

		# test size
		raise "Wrong size of vector" unless vector.size == size * 3
		# copy values
		each_index {|i|
			atom = at(i)
			atom.x = vector[i*3]
			atom.y = vector[i*3+1]
			atom.z = vector[i*3+2]
		}
		return nil
	end

	#=======================================================================
	# Global properties
	#=======================================================================
	
	def mass
		m = 0.0
		each {|atom| m += atom.mass}
		return m
	end
	
	def moment_of_inertia # => Matrix
		# move center of mass to [0,0,0]
		com = center_of_mass

		# 3x3 array
		i = Matrix.zero(3)

		each {|atom|
			i[0,0] += atom.mass * ( (atom.y-com.y)**2 + (atom.z-com.z)**2)		# xx
			i[1,1] += atom.mass * ( (atom.x-com.x)**2 + (atom.z-com.z)**2)		# yy
			i[2,2] += atom.mass * ( (atom.x-com.x)**2 + (atom.y-com.y)**2)		# zz

			i[1,0] -= atom.mass * (atom.x-com.x) * (atom.y-com.y)
			i[2,0] -= atom.mass * (atom.x-com.x) * (atom.z-com.z)
			i[2,1] -= atom.mass * (atom.y-com.y) * (atom.z-com.z)
		}
		i[0,1] = i[1,0]
		i[0,2] = i[2,0]
		i[1,2] = i[2,1]

		return i
	end

	def moment_of_inertia_diagonalized # => [eigenvectors:[Vectors], eigenvalues:[numbers]]
		mmt = moment_of_inertia
		vec, real , imag = mmt.eigensystem
		Matrix.sort_eigensystem!(vec, real , imag, :ascending_real)

		real = real.to_a.flatten # eigenvalues as array
		vectors = [] # eigenvectors as array of vectors
		3.times{|i| vectors[i] = vec.column_as_vector(i)}

		# Handedness fix - make sure the three molecular axes are right-handed
		v1, v2, v3 = vectors
		if (v1.cross_product(v2)).dot(v3) < 0
			vectors[2] *= -1.0
		end

		return [vectors,real]
	end

	#--------------------------------------------------
	# def linear?(threshold = 0.01) # => Boolean
	# 	#: Determines whether the geometry is linear, within the given tolerance (in A)
	# 	# Diatomics are always linear
	# 	return true if size <= 2
	# 	# Find two most distant atoms
	# 	a1 = nil
	# 	a2 = nil
	# 	maxdistance = 0.0
	# 	each_index{|i|
	# 		i.times{|j|
	# 			d = at(i).distance(at(j))
	# 			if d > maxdistance
	# 				maxdistance = d
	# 				a1 = at(i)
	# 				a2 = at(j)
	# 			end
	# 		}
	# 	}
	# 	# Measure distance of other atoms from the axis
	# 	each_with_index{|atom,i|
	# 		next if a1 == atom || a2 == atom
	# 		return false if atom.point_line_distance(a1,a2) > threshold
	# 	}
	# 	# All atoms are on the axis
	# 	return true
	# end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# def linear2? #=> Boolean
	# 	#: Determines whether the geometry is linear (with some numerical tolerance).
	# 	#: Should be faster than "linear?"
	# 	moment_of_inertia_diagonalized[1].each{|eigenval|
	# 		return true if eigenval.abs < 1.0e-7
	# 	}
	# 	return false
	# end
	#-------------------------------------------------- 

end
