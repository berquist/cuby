################################################################################
#
# class Vector
#
# Author: Jan Rezac
# Date created: 2008-08-30
# License: Cuby license
# Description: Vector algebra
# Status: Tested
#
################################################################################

#! <h1>class Vector</h1>

#: Cuby's implementation of vector algebra and basic 3D geometry
#: in Vector class.

#: This is pure ruby implementation, in contrast to math_vector which uses
#: binary library. It is intended as a replacement where the performance
#: is not crucial, but portability is.

module Algebra
class Vector
	@@default_epsilon = 1.0e-8

	#=======================================================================
	# initializators
	#=======================================================================

	def self.[](*array) # new instance of the class
		#: Creates new vector using array notation (without brackets)
		#| v = Vector[1.0, 0.0, 0.0]
		return from_array(array)
	end

	def self.from_values(*array) # new instance of the class
		#: Creates new vector from list of values
		#| v = Vector.from_values(1.0, 2.0, 3.0)
		return from_array(array)
	end

	def self.random(size)
		#: Creates a new vector filled with random values 0.0 <= r < 1.0
		x = self.of_size(size)
		x.each_index{|i|
			x[i] = rand
		}
		return x
	end

	#=======================================================================
	# Class functions
	#=======================================================================
	
	def self.orthonormalization_gram_schmidt(array_of_vectors) # => nil
		#: Trannsforms set of vector into orthogonal set of normalized vectors
		#: using the Gram-Schmidt algorithm
		array_of_vectors.size.times {|j|
			j.times{|i|
				array_of_vectors[j] = array_of_vectors[j] -
					((array_of_vectors[i].dot(array_of_vectors[j])) /
					 (array_of_vectors[i].dot(array_of_vectors[i]))) *
					 array_of_vectors[i]
			}
			array_of_vectors[j].normalize!
		}
		return nil
	end

	#=======================================================================
	# Access
	#=======================================================================
	
	def first # => Float
		#: Return first element of the vector
		return self[0]
	end

	def last # => Float
		#: Return last element of the vector
		return self[size - 1]
	end
	
	#=======================================================================
	# basics
	#=======================================================================

	def copy_from!(source) # nil
		#: Copies value from source Vector(or derivative) or Array, conserving the object id.
		#| v = Vector[0,1,2]
		#| v.copy_from [4,5,6]
		#| => Vector: [4.0, 5.0, 6.0]
		
		raise(DimensionError,"Copy_from can copy only from objects with the same size.") if size != source.size
		source.each_index {|i|
			item = source[i]
			raise(TypeError,"Vector can be created only from array of numbers") unless item.kind_of?(Numeric)
			self[i] = item
		}
		return nil
	end

	def clone # Vector
		#: Synonym for -->deep_copy
		return self.deep_copy
	end

	#=======================================================================
	# type conversion
	#=======================================================================
	
	def to_s # String
		#: Conversion to string
		#| Vector[1.0, 0.0, 0.0].to_s
		#| => "[1.0, 0.0, 0.0]"
		return self.class.to_s + ": [" + self.to_a.flatten.join(", ") + "]"
		#return @vector_data_struct.to_a.flatten.join("\n")
	end

	def inspect # String
		#: Equal to -->to_s with prepended newline for nice irb output.
		return "\n" + to_s
	end

	def to_vector # Vector
		#: Allows any descendants to be converted to vector
		return self.class.from_array(self.to_a)
	end

	#=======================================================================
	# iterators
	#=======================================================================
	
	def each_index
		#: Iterate over indexes
		size.times {|i| yield(i)}
	end

	def each
		#: Iterate over elements
		size.times {|i| 
			yield self[i]
		}
	end

	def each_with_index
		#: Iterate over elements, yielding both value and index
		size.times {|i| yield(self[i], i)}
	end
	
	def collect # Array
		#: Iterate over elements, results of the block are returned as an array
		result = []
		size.times {|i| 
			result[i] = yield self[i]
		}
		return result
	end

	#=======================================================================
	# unary operators
	#=======================================================================
	
	def +@ # self
		#: Unary plus. Does nothing.
		return self
	end

	def -@ # instance of self.class
		#: Unary minus. Equal to self * -1.0
		return self * -1.0
	end
	
	#=======================================================================
	# operators: vector - vector
	#=======================================================================

	def dot_product(vector) # Float
		#: Dot product
		dot(vector)
	end

	#=======================================================================
	# operators: vector - something else
	#=======================================================================
	
	def *(operand) # instance of self.class | Float (Vector * Vector)
		#: Multiplication operator, allows vector * Float, Fixnum or Vector
		if operand.kind_of?(Numeric)
			return self.comutative_multiply(operand)
		elsif operand.kind_of?(Vector)
			return self.dot(operand)
		else
			raise(TypeError, "Vector.*: Can't multiply vector by class #{operand.class}")
		end
	end

	def /(operand) # instance of self.class
		#: Division of vector by numbers
		raise(TypeError,"Vector can be divided only by a number") unless operand.kind_of?(Numeric)
		return self * (1.0 / operand)
	end

	#=======================================================================
	# comparison
	#=======================================================================
	
	def self.epsilon
		#: Returns threshold for loose comparison using -->=~ operator.
		return @@default_epsilon
	end

	def self.epsilon=(value)
		#: Sets threshold for loose comparison using -->=~ operator.
		@@default_epsilon = value
	end

	#=======================================================================
	# misc
	#=======================================================================
	
	def absolute # Float
		#: Synonym to -->abs
		abs
	end

	def r # Float
		#: Synonym to -->abs
		abs
	end

	def distance(vector) # Float
		#: Distance between two points specified by vectors
		return (self - vector).absolute
	end

	def dist(vector) # Float
		#: Synonym for -->distance
		distance(vector)
	end

	def rms
		return abs/size**0.5
	end

	#=======================================================================
	# YAML support
	#=======================================================================
	
	require "yaml"

	begin
		yaml_engine = "psych" if YAML == Psych
		yaml_engine = "syck" if YAML == Syck
	rescue
		# ruby1.8 uses syck by default
		yaml_engine = 'syck'
	end

	if yaml_engine == 'syck'
		# Syck yaml engine
	
		yaml_as "tag:cuby.molecular.cz,2009:#{self}"

		def to_yaml(opts = {})
			YAML::quick_emit( self.object_id, opts ) do |out|
				out.map(taguri) do |map|
					a = self.to_a
					def a.to_yaml_style; :inline; end # make it inline
					map.add('elements', a)
				end
			end
		end

		def self.yaml_new(klass, tag, val)
			array = val['elements']
			return Vector.from_array(array)
		end

	elsif yaml_engine == 'psych'
		# Psych yaml engine
	
		def encode_with(coder)
			coder.tag = "!cuby.molecular.cz,2009/Algebra::Vector"
			arr = self.to_a
			coder.style = Psych::Nodes::Mapping::FLOW # BLOCK or FLOW
			coder['elements'] = arr
		end

		def init_with(coder)
			puts "init_with"
		end

		YAML::add_domain_type('cuby.molecular.cz,2009', 'Algebra::Vector') { |type, value|
			Vector.from_array(value['elements'])
		}
	else
		raise "YAML engine can not be determined"
	end
	
end

#!
class Not3DError < RuntimeError
end
#!!

module Vector3DMethods

#: Additional functionality available only in 3D vectors
#: Module Vector3D is automatically included in Vector class

	# check size of the vector, raise exception if size != 3
	def raise_if_not_3D
		#: Internal use, check for size == 3, raises exception otherwise
		raise(Not3DError,"This operation is available only for 3D vectors.") if size != 3
	end

	def cross_product(vector) # Vector
		# Cross product, implemented for 3D vectors only
		raise_if_not_3D
		vector.raise_if_not_3D
		
		x = self[1] * vector[2] - self[2] * vector[1]
		y = self[2] * vector[0] - self[0] * vector[2]
		z = self[0] * vector[1] - self[1] * vector[0]

		return self.class[x,y,z]
	end

	def angle(vector) # Float
		#: Angle between two 3D vectors (in radians)
		raise_if_not_3D
		return Math.acos(self.dot(vector)/self.absolute/vector.absolute)
	end

	def angle_deg(vector) # Float
		#: Angle between two 3D vectors (in degrees)
		return self.angle(vector) / Math::PI * 180.0
	end

	def rotate_by_matrix(rotator) # Vector
		#: Rotation of vector by application of 3x3 matrix
		raise_if_not_3D

		v2 = (rotator*self)
		return v2.to_vector
	end

	def rotate_by_matrix!(rotator) # nil
		#: In-place rotation of vector by application of 3x3 matrix
		raise_if_not_3D

		v2 = (rotator*self)
		self.copy_from! v2.to_vector
		return nil
	end

	def self.rotator_from_angle_axis(angle,axis) # Matrix
		#: Creates rotator matrix from ange (float) and axis (3D vector)

		axis.raise_if_not_3D

		axis = axis.normalize

		x = axis[0]
		y = axis[1]
		z = axis[2]

		cosa = Math.cos(angle) 
		sina = Math.sin(angle) 

		r = Matrix[	[cosa + (1-cosa)*x*x,	(1-cosa)*x*y - sina*z,	(1-cosa)*x*z + sina*y],
				[(1-cosa)*x*y + sina*z,	cosa + (1-cosa)*y*y,	(1-cosa)*z*y - sina*x],
				[(1-cosa)*x*z - sina*y,	(1-cosa)*z*y + sina*x,	cosa + (1-cosa)*z*z]]

		return r
	end

	def rotate_angle_axis(angle,axis) # Vector
		#: Returns vector rotated by angle (float) around axis (3D vector)
		raise_if_not_3D
		axis.raise_if_not_3D

		r = Vector3DMethods.rotator_from_angle_axis(angle,axis)
		return rotate_by_matrix(r)
	end

	def rotate_angle_axis!(angle,axis) # nil
		#: Rotates vector by angle (float) around axis (3D vector)
		raise_if_not_3D
		axis.raise_if_not_3D

		r = Vector3DMethods.rotator_from_angle_axis(angle,axis)
		rotate_by_matrix!(r)
		return nil
	end

	#=======================================================================
	# Direct access to coordinates
	#=======================================================================
	
	def x # Float
		raise_if_not_3D
		self[0]
	end
	def y # Float
		raise_if_not_3D
		self[1]
	end
	def z # Float
		raise_if_not_3D
		self[2]
	end

	def x=(value)
		raise_if_not_3D
		self[0] = value
		return value
	end

	def y=(value)
		raise_if_not_3D
		self[1] = value
		return value
	end

	def z=(value)
		raise_if_not_3D
		self[2] = value
		return value
	end
end

#!
class Vector
	include Vector3DMethods
end
#!!
end
