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

#: Access to the object can be locked using Locking module,
#: all methods are locked

#: This is pure ruby implementation, in contrast to math_vector which uses
#: binary library. It is intended as a replacement where the performance
#: is not crucial, but portability is.

module Algebra

#!
class DimensionError < RuntimeError
end
##!

class Vector
	#internal data storage: @vector_data_struct = Array

	#!
	def initialize(array)
		# This method should not be called directly. Internally, the argument should be an
		# array with elements of the vector. However, for compatibility with the C version,
		# it accepts an integer and it builds a zero vector of that size.
		unless array.class == Array
			if array.kind_of?(Integer)
				array = Array.new(array, 0.0)
			else
				raise(TypeError,"Vector.new argument must be array or integer.
(In general, Vector.new should not be used outside the Vector class,
use Vector.zero or Vector.from_array instead.)")
			end
		end
		@vector_data_struct = array
	end

	def get_data_struct # Linalg::DMatrix
		# Internal use only! Do not use outside this file 
		return @vector_data_struct
	end

	private :get_data_struct

	#!!
	
	#=======================================================================
	# initializators
	#=======================================================================

	def self.from_array(array) # new instance of the class
		#: Creates new vector from array
		#| v1 = Vector.from_array([1.0, 0.0, 0.0])
		#| v2 = Vector.from_array(some_array)
		raise(TypeError,"Vector can not be created from an empty array") if array.size == 0
		a = []
		array.each {|item|
			raise(TypeError,"Vector can be created only from array of numbers") unless item.kind_of?(Numeric)
			a << item.to_f
		}
		return self.new(a)
	end

	def self.of_size(size,fill=0.0) # new instance of the class
		#: Creates new vector of given size filled with number fill
		#| Vector.of_size(5,1.0)
		#| => Vector: [1.0,1.0,1.0,1.0,1.0]
		raise(TypeError,"Vector can not have size 0 or smaller") if size <= 0
		return self.new(Array.new(size,fill))
	end

	def self.zero(size) # new instance of the class
		#: Creates new vector of given size filled with zeroes
		#| Vector.of_size(5)
		#| => Vector: [0.0,0.0,0.0,0.0,0.0]
		raise(TypeError,"Vector can not have size 0 or smaller") if size <= 0
		return self.new(Array.new(size,0.0))
	end

	#=======================================================================
	# basics
	#=======================================================================
	
	def size # Fixnum
		#: Returns size of the vector
		#| Vector[1.0, 0.0, 0.0].size
		#| => 3
		return @vector_data_struct.size
	end

        def deep_copy # Vector
		#: deep copy
		d = @vector_data_struct.dup
		return self.class.new(d)
	end

	#=======================================================================
	# element access
	#=======================================================================
	
	def [](index) # Float
		#: Returns selected element
		raise(IndexError, "Index out of bounds") if index >= @vector_data_struct.size || index < 0
		return @vector_data_struct[index]
	end

	def []=(index,value) # Float (value written)
		#: Write to selected element
		raise(IndexError, "Index out of bounds") if index >= @vector_data_struct.size || index < 0
		raise(TypeError) unless value.kind_of?(Numeric)
		@vector_data_struct[index] = value.to_f
		return value
	end

	def subvector(first, size) # => Vector
		v = self.class.zero(size)
		size.times {|i|
			v[i] = self[first + i]
		}
		return v
	end
	
	#=======================================================================
	# type conversion
	#=======================================================================
	
	def to_a # Array of Floats
		#: Converts Vector to Array
		return @vector_data_struct.dup
	end

	def to_matrix # Matrix
		#: Converts Vector to column Matrix
		return Matrix.columns([self.to_a])
	end


	#=======================================================================
	# iterators
	#=======================================================================
	
	#=======================================================================
	# unary operators
	#=======================================================================
	
	#=======================================================================
	# operators: vector - vector
	#=======================================================================
	
	def +(vector) # Vector
		#: Vector addition
		raise(TypeError, "Argument is not a Vector") unless vector.kind_of?(Vector)
		# size check
		raise(DimensionError,"Can't add vectors of different sizes") if size != vector.size
		result = []
		@vector_data_struct.each_index{|i|
			result[i] = @vector_data_struct[i] + vector[i]
		}
		return self.class.new(result)
	end

	def -(vector) # Vector
		#: Vector subtraction
		raise(TypeError, "Argument is not a Vector") unless vector.kind_of?(Vector)
		# size check
		raise(DimensionError,"Can't subtract vectors of different sizes") if size != vector.size
		result = []
		@vector_data_struct.each_index{|i|
			result[i] = @vector_data_struct[i] - vector[i]
		}
		return self.class.new(result)
	end

	def plus!(vector) # nil
		#: In-place vector addition, faster than using a = a + b
		raise(TypeError, "Argument is not a Vector") unless vector.kind_of?(Vector)
		# size check
		raise(DimensionError,"Can't add vectors of different sizes") if size != vector.size
		@vector_data_struct.each_index{|i|
			@vector_data_struct[i] = @vector_data_struct[i] + vector[i]
		}
		return nil
	end

	def minus!(vector) # nil
		#: In-place vector subtraction, faster than using a = a - b
		raise(TypeError, "Argument is not a Vector") unless vector.kind_of?(Vector)
		# size check
		raise(DimensionError,"Can't subtract vectors of different sizes") if size != vector.size
		result = []
		@vector_data_struct.each_index{|i|
			@vector_data_struct[i] = @vector_data_struct[i] - vector[i]
		}
		return nil
	end

	def dot(vector) # Float
		#: Dot product, synonym to -->dot_product
		raise(TypeError, "Argument is not a Vector") unless vector.kind_of?(Vector)
		# size check
		raise(DimensionError,"Can't multiply vectors of different sizes") if size != vector.size
		result = 0.0
		@vector_data_struct.each_index{|i|
			result += @vector_data_struct[i] * vector[i]
		}
		return result
	end

	def elementwise_multiply(vector) # => Vector
		raise(TypeError, "Argument is not a Vector") unless vector.kind_of?(Vector)
		# size check
		raise(DimensionError,"Can't multiply vectors of different sizes") if size != vector.size
		result = self.class.zero(vector.size)
		@vector_data_struct.each_index{|i|
			result[i] = @vector_data_struct[i] * vector[i]
		}
		return result
	end

	def elementwise_divide(vector) # => Vector
		raise(TypeError, "Argument is not a Vector") unless vector.kind_of?(Vector)
		# size check
		raise(DimensionError,"Can't divide vectors of different sizes") if size != vector.size
		result = self.class.zero(vector.size)
		@vector_data_struct.each_index{|i|
			result[i] = @vector_data_struct[i] / vector[i]
		}
		return result
	end

	#=======================================================================
	# operators: vector - something else
	#=======================================================================
	
	def comutative_multiply(operand)
		#: this method is equivalent to * operator, but is used by the
		#: numerical types extension found in algebra.rb
		result = self.class.zero(size)
		@vector_data_struct.each_index{|i|
			result[i] = @vector_data_struct[i] * operand
		}
		return result
	end

	#=======================================================================
	# comparison
	#=======================================================================
	
	def == (vector) # true | false
		# Comparison operator
		if vector.kind_of?(Vector)
			raise(DimensionError,"Vectors of different sizes can't be compared") unless size == vector.size
			@vector_data_struct.each_index{|i| return false if @vector_data_struct[i] != vector[i]}
			return true
		else
			return false
		end
	end

	def =~ (vector) # true | false
		#: Loose comparison of two vectors, element by element. Max. allowed
		#: error is stored in class variable, accessible by methods
		#: -->self.epsilon and -->self.epsilon= .
		if vector.kind_of?(Vector)
			raise(DimensionError,"Vectors of different sizes can't be compared") unless size == vector.size
			@vector_data_struct.each_index{|i| return false if (@vector_data_struct[i] - vector[i]).abs > @@default_epsilon}
			return true
		else
			return false
		end
	end

	#=======================================================================
	# misc
	#=======================================================================
	
	def sum
		sum = 0
		each {|e| sum += e}
		return sum
	end

	def abs # Float
		#: Absolute value (size) of vector
		sum = 0.0
		each {|e| sum += e**2}
		return sum**0.5
	end

	def max
		max = self[0]
		each {|e| max = e if e > max}
		return max
	end

	def min
		min = self[0]
		each {|e| min = e if e < min}
		return min
	end

	def max_abs
		max = self[0].abs
		each {|e| max = e.abs if e.abs > max}
		return max
	end

	def min_abs
		min = self[0].abs
		each {|e| min = e.abs if e.abs < min}
		return min
	end

	def normalize # Vector
		#: Returns normalized vector
		a = abs
		return self / a
	end

	def normalize! # nil
		#: In-place version of -->normalize
		a = abs
		@vector_data_struct.each_index {|i| @vector_data_struct[i] = @vector_data_struct[i] / a}
		return nil
	end

end
end
