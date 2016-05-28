################################################################################
#
# class Coordinate - basics
#
# Author: Jan Rezac
# Date created: 2008-08-30
# License: Cuby license
# Description: Coordinate in 3D space
# Status: Documented, Tested
#
################################################################################

#: Basic functionality of Coordinate class. Can be replaced by C extension.
#: This file should not be loaded directly, but via coordinate.rb

class Coordinate
	@@epsilon = 1.0e-8

	attr_reader :x
	attr_reader :y
	attr_reader :z

	def Coordinate.version
		return :ruby
	end

	#> Coordinate.new(x = 0.0, y = 0.0, z = 0.0) # => Coordinate
	#!
	def initialize(x = 0.0, y = 0.0, z = 0.0)
	#!!
		#: Creates new Coordinate from three values, ommited arguments
		#: are set to zero
		if x.class != Float
			if x.kind_of? Numeric
				x = x.to_f
			else
				raise(TypeError, "Coordinate element must be a number (is #{x.class})")
			end
		end
		if y.class != Float
			if y.kind_of? Numeric
				y = y.to_f
			else
				raise(TypeError, "Coordinate element must be a number (is #{y.class})")
			end
		end
		if z.class != Float
			if z.kind_of? Numeric
				z = z.to_f
			else
				raise(TypeError, "Coordinate element must be a number (is #{z.class})")
			end
		end
		@x = x
		@y = y
		@z = z
	end
	
	#=======================================================================
	# initializators
	#=======================================================================

	def self.[](*array) # => Coordinate
		#: Creates new Coordinate using array notation (without brackets)
		#| c = Coordinate[1.0, 0.0, 0.0]
		return Coordinate.from_array(array)
	end

	def self.from_array(array) # => Coordinate
		#: Creates new coordinate from array of three floats
		#| c1 = Coordinate([1.0, 0.0, 0.0])
		#| c2 = Coordinate(some_array)
		raise(TypeError,"Coordinate can be created only from 3-dimensional array") if array.size != 3
		return Coordinate.new(array[0],array[1],array[2])
	end

	def self.from_coordinate(coordinate) # => Coordinate
		# Creates coordinate from coordinate or its descendant
		raise(TypeError,"Coordinate can be created only from a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		return Coordinate.new(coordinate.x, coordinate.y, coordinate.z)
	end

	#=======================================================================
	# Write to components
	#=======================================================================

	#: Following three methods provide write to the coordinate elements
	#: with type checking.
		
	def x=(value) # => Float (value written)
		raise(TypeError,"Coordinate element must be Float") if value.class != Float 
		@x = value
	end
	
	def y=(value) # => Float (value written)
		raise(TypeError,"Coordinate element must be Float") if value.class != Float 
		@y = value
	end
	
	def z=(value) # => Float (value written)
		raise(TypeError,"Coordinate element must be Float") if value.class != Float 
		@z = value
	end

	#=======================================================================
	# basics
	#=======================================================================

	def deep_copy # => Coordinate
		#: Creates real copy of the object
		return Coordinate.new(@x, @y, @z)
	end
	
	def dup # => Coordinate
		#: Creates real copy of the object
		return Coordinate.new(@x, @y, @z)
	end

	def set_coord(coordinate) # => nil
		#: Sets the coordinate (its elements) from other object
		raise(TypeError,"Coordinates can be set ofnly from a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		self.x = coordinate.x
		self.y = coordinate.y
		self.z = coordinate.z
		return nil
	end

	#=======================================================================
	# element access
	#=======================================================================
	
	def [](index) # => Float
		#: Returns selected element
		#| c = Coordinate[1.0, 2.0, 3.0]
		#| c[2]
		#| => 3.0
		case index
		when 0
			return @x
		when 1
			return @y
		when 2
			return @z
		else
			raise(IndexError,"Index not in range 0..2")
		end
	end

	def []=(index,value) # => Float (value written)
		#: Write to element selected by index
		#| c = Coordinate[0.0, 0.0, 0.0]
		#| c[0] = 1.0
		#| c
		#| => Coordinate: [1.0, 0.0, 0.0]
		case index
		when 0
			self.x = value
		when 1
			self.y = value
		when 2
			self.z = value
		else
			raise(IndexError,"Index not in range 0..2")
		end
		return value
	end

	#=======================================================================
	# type conversion
	#=======================================================================
	
	def to_s # => String
		#: Conversion to string
		#| Coordinate[1.0, 2.0, 3.0].to_s
		#| => "Coordinate: [1.0, 2.0, 3.0]"
		return self.class.to_s + ": [#{@x}, #{@y}, #{@z}]"
	end

	def inspect # => String
		#: Equal to -->to_s with prepended newline for nice irb output.
		return "\n" + to_s
	end

	def to_a # => Array of Floats
		#: Converts Coordinate to Array
		#| Coordinate[1.0, 2.0, 3.0].to_a
		#| => [1.0, 2.0, 3.0]
		return [@x,@y,@z]
	end

	def to_coordinate # => Coordinate
		#: Allows conversion of descendants to Coordinate
		return Coordinate.new(@x, @y, @z)
	end

	#=======================================================================
	# iterators
	#=======================================================================
	
	def each_index # => nil
		#: Iterate over indexes
		3.times {|i| yield(i)}
		return nil
	end

	def each # => nil
		#: Iterate over elements
		yield(@x)
		yield(@y)
		yield(@z)
		return nil
	end
	
	def collect # => Coordinate
		#: Iterate over elements, results of the block are returned as new Coordinate
		xx = yield(@x)
		yy = yield(@y)
		zz = yield(@z)
		return Coordinate.new(xx, yy, zz)
	end

	def map # => Coordinate
		#: Synonym for -->collect
		return collect
	end

	def map! # => nil
		#: Iterate over elements, result of the block replaces original value of element
		self.x= yield(@x)
		self.y= yield(@y)
		self.z= yield(@z)
		return nil
	end

	#=======================================================================
	# unary operators
	#=======================================================================
	
	def +@ # => self
		#: Unary plus. Does nothing.
		return self
	end

	def -@ # => Coordinate
		#: Unary minus. Returns deep copy with all x,y,z multiplied by -1.0
		result = self.deep_copy
		result.x *= -1.0
		result.y *= -1.0
		result.z *= -1.0
		return result
	end
	
	#=======================================================================
	# operators: coordinate - coordinate
	#=======================================================================
	
	def +(coordinate) # => Coordinate
		#: Coordinate addition
		raise(TypeError, "Operand must be a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		return Coordinate.new(@x+coordinate.x, @y+coordinate.y, @z+coordinate.z)
	end

	def -(coordinate) # => Coordinate
		#: Coordinate subtraction
		raise(TypeError, "Operand must be a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		return Coordinate.new(@x-coordinate.x, @y-coordinate.y, @z-coordinate.z)
	end

	def plus!(coordinate) # => nil
		#: In-place coordinate addition, faster than using a = a + b
		raise(TypeError, "Operand must be a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		@x += coordinate.x
		@y += coordinate.y
		@z += coordinate.z
		return nil
	end

	def minus!(coordinate) # => nil
		#: In-place coordinate subtraction, faster than using a = a - b
		raise(TypeError, "Operand must be a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		@x -= coordinate.x
		@y -= coordinate.y
		@z -= coordinate.z
		return nil
	end

	def dot_product(coordinate) # => Float
		#: Dot product
		raise(TypeError, "Operand must be a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		return @x * coordinate.x + @y * coordinate.y + @z * coordinate.z
	end

	def dot(coordinate) # => Float
		#: Dot product, synonym for -->dot_product
		dot_product(coordinate)
	end

	def cross_product(coordinate) # => Coordinate
		#: Cross product
		raise(TypeError, "Operand must be a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		rx = @y*coordinate.z - @z*coordinate.y
		ry = @z*coordinate.x - @x*coordinate.z
		rz = @x*coordinate.y - @y*coordinate.x
		return Coordinate.new(rx,ry,rz)
	end

	#=======================================================================
	# operators: Coordinate - something else
	#=======================================================================
	
	def *(operand) # => Coordinate | Float (Coordinate * Coordinate)
		#: Multiplication operator, allows Coordinate * Numeric or Coordinate,
		#: returns dot product for Coordinate * Coordinate.
		if operand.kind_of?(Numeric)
			return Coordinate.new(@x*operand,@y*operand,@z*operand)
		elsif operand.kind_of?(Coordinate)
			return dot_product(operand)
		else
			raise(TypeError, "Can't multiply Coordinate by class #{operand.class}")
		end
	end

	def /(operand) # => Coordinate
		#: Division of coordinate by number (Numeric)
		if operand.kind_of?(Numeric)
			return Coordinate.new(@x/operand,@y/operand,@z/operand)
		else
			raise(TypeError, "Can't divide Coordinate by class #{operand.class}")
		end
	end

	#=======================================================================
	# comparison
	#=======================================================================
	
	def == (coordinate) # => true | false
		# Comparison operator
		return false unless coordinate.kind_of?(Coordinate)
		return (@x==coordinate.x && @y==coordinate.y && @z==coordinate.z)
	end

	def =~ (coordinate) # => true | false
		#: Loose comparison, difference between each pair of elements must be
		#: smaller than class variable epsilon.
		return false unless coordinate.kind_of?(Coordinate)
		return ((@x-coordinate.x).abs < @@epsilon && 
		        (@y-coordinate.y).abs < @@epsilon &&
		        (@z-coordinate.z).abs < @@epsilon )
	end

	def zero?
		# Return true if the coordinate is exactly [0,0,0]
		return @x == 0.0 && @y == 0.0 && @z == 0.0
	end

	def Coordinate.epsilon=(value) # => Float (value written)
		#: Sets threshold for loose comparison (class variable epsilon)
		@@epsilon = value
	end

	def Coordinate.epsilon # => Float
		#: Reads threshold for loose comparison (class variable epsilon)
		return @@epsilon
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def absolute # => Float
		# Absolute value (size) of vector
		return (@x**2 + @y**2 + @z**2)**0.5
	end

	def abs # => Float
		#: Synonym for -->absolute
		absolute
	end

	def r # => Float
		#: Synonym for -->absolute
		absolute
	end

	def normalize # => Coordinate
		#: Returns normalized vector
		a = absolute
		return Coordinate.new(@x/a,@y/a,@z/a)
	end

	def normalize! # => nil
		#: In-place version of -->normalize
		a = absolute
		@x /= a
		@y /= a
		@z /= a
		return nil
	end

	def distance(coordinate) # => Float
		#: Distance between two points specified by coordinates
		raise(TypeError, "Argument must be a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		return (self - coordinate).absolute
	end

	def dist(coordinate) # => Float
		#: Synonym for -->distance
		distance(coordinate)
	end

	def angle(coordinate) # => Float
		#: Angle between two 3D vectors in space (in radians)
		raise(TypeError, "Argument must be a Coordinate-like type") unless coordinate.kind_of?(Coordinate)
		raise(ArgumentError, "Coordinate for angle calculation can not be zero") if self.zero?
		raise(ArgumentError, "Coordinate for angle calculation can not be zero") if coordinate.zero?
		cos = self.dot_product(coordinate)/self.absolute/coordinate.absolute
		# Handle numerical issues
		cos = -1.0 if cos < -1.0
		cos = 1.0 if cos > 1.0
		return Math::acos(cos)
	end

	def angle_deg(coordinate) # => Float
		#: Angle between two 3D vectors in space (in degrees)
		return self.angle(coordinate) / Math::PI * 180.0
	end

end

