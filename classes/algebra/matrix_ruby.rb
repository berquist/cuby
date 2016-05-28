################################################################################
#
#  class Matrix
#
################################################################################

#! <h1>class Matrix</h1>

#: Cuby's linear algebra library.
#: This is pure ruby implementation - computationaly intensive operations won't
#: be efficient.
#: If efficiency is required, use the binary extension

module Algebra
require @@directory+"/matrix_ruby/matrix_ruby_eigensystem"
require @@directory+"/matrix_ruby/matrix_ruby_LU"
require @@directory+"/matrix_ruby/matrix_ruby_QR"
require @@directory+"/matrix_ruby/matrix_ruby_SVD"

# Exceptions
class  DimensionError < RuntimeError
end

class Matrix

	#!
	#=======================================================================
	# Internal use only
	#=======================================================================
	
	def initialize(data, second=nil)
		# creates new matrix, taking an array of arrays as the internal matrix representation
		# the data are not copied!
		# there is no error check!
		# this must be provided in the method calling .new
	
		# This method should not be called directly. Internally, the argument should be an
		# array of arrays containing the elements. However, for compatibility with the C version,
		# it accepts also two integers and it builds a zero matrix of that size
		unless data.class == Array
			if data.kind_of?(Integer) && second.kind_of?(Integer)
				m = data
				n = second
				data = []
				m.times {
					data << Array.new(n,0.0)
				}
			else
				raise(TypeError,"Matrix.new argument(s) must be either
one array or two integers.
(In general, Matrix.new should not be used outside the Matrix class,
use Matrix.zero or Matrix.from_array instead.)")
			end
		end
		@matrix_data_struct = data
	end

	def get_data_struct
		return @matrix_data_struct
	end
	#!!

	#=======================================================================
	# initializators
	#=======================================================================

	def self.from_array(array)
		#: Creates matrix from array of arrays (used with brackets)
		#| a = Matrix([[1,2],[3,4]])
		#| =>
		#| |1.0 2.0|
		#| |3.0 4.0|
		raise(TypeError,"Matrix can be created only from Array") unless array.class == Array
		raise(TypeError,"Matrix can be created only from non-empty Array") unless array.size > 0
		raise(TypeError,"Matrix can be created only from Array of Arrays") unless array[0].class == Array
		raise(TypeError,"Matrix can be created only from non-empty Arrays") unless array[0].size > 0
		a = []
		firstrowsize = array[0].size
		array.each {|row| 
			raise(TypeError,"Matrix can be created only from Arrays of the same size") unless row.size == firstrowsize
			newrow = []
			row.each {|item|
				raise(TypeError,"Matrix can be created only from 2D Array of numbers") unless item.kind_of?(Numeric)
				newrow << item.to_f
			}
			a << newrow
		}
		return self.new(a)
	end

	def self.filled(m,n, fill)
		#: Creates zero-filled matrix
		raise(TypeError,"Matrix dimension m must be larger than 0") if m <= 0
		raise(TypeError,"Matrix dimension n must be larger than 0") if n <= 0
		data = []
		m.times {
			data << Array.new(n,fill)
		}
		return self.new(data)
	end

	def self.zero(m,n=nil)
		#: Creates zero-filled matrix
		#: if n is ommitted, square matrix is created
		n = m if n == nil
		return self.filled(m,n, 0.0)
	end

	def self.diagonal1(size, value)
		newmatrix = self.zero(size)
		newmatrix.n.times{|i|
			newmatrix[i,i] = value
		}
		return newmatrix
	end

	def self.diagonal2(array)
		newmatrix = self.zero(array.size)
		newmatrix.n.times{|i|
			raise(TypeError,"Matrix can be constructed only from Array of numbers") unless array[i].kind_of?(Numeric)
			newmatrix[i,i] = array[i].to_f
		}
		return newmatrix
	end

	#=======================================================================
	# Version
	#=======================================================================
	
	def self.version
		return :ruby
	end

	#=======================================================================
	# basics
	#=======================================================================
	
	def m
		#: Vertical size (number of rows)
		return @matrix_data_struct.size
	end

	def n
		#: Horizontal size (number of columns)
		return @matrix_data_struct[0].size
	end

	def deep_copy
		#: Deep copy
		dataclone = []
		@matrix_data_struct.each {|row| dataclone << row.dup}
		return self.class.new(dataclone)
	end

	#=======================================================================
	# element access
	#=======================================================================
	
	def [](row,col)
		#: Returns selected element
		raise(IndexError,"Row index (i) out of bounds") if row >= m || row < 0
		raise(IndexError,"Col index (j) out of bounds") if col >= n || col < 0
		return @matrix_data_struct[row][col]
	end

	def []=(row,col,value)
		#: Write to selected element
		raise(IndexError,"Row index (i) out of bounds") if row >= m || row < 0
		raise(IndexError,"Col index (j) out of bounds") if col >= n || col < 0
		raise(TypeError,"Only numbers can be stored in matrix") unless value.kind_of?(Numeric)
		@matrix_data_struct[row][col] = value.to_f
	end

	def submatrix(i, j, size_m, size_n)
		#: Returns selected region of the matrix
		new = self.class.zero(size_m, size_n)
		size_m.times {|ii|
			size_n.times {|jj|
				new[ii,jj] = self[ii+i,jj+j]
			}
		}
		return new
	end

	def row_as_array(index)
		#: Returns selected row as an array
		a = []
		n.times {|j|
			a << self[index,j]
		}
		return a
	end

	def column_as_array(index)
		#: Returns selected column as an array
		a = []
		m.times {|j|
			a << self[j,index]
		}
		return a
	end

	def col_as_array(index)
		return column_as_array(index)
	end

	def row_as_vector(index)
		#: Returns selected row as a Vector
		return Vector.from_array(row_as_array(index))
	end

	def column_as_vector(index)
		#: Returns selected column as a Vector
		return Vector.from_array(column_as_array(index))
	end

	def col_as_vector(index)
		return column_as_vector(index)
	end

	def diagonal_as_vector
		# Returns the diagonal of the matrix as a vector
		dia_size = [m,n].min
		vec = Vector.zero(dia_size)
		dia_size.times {|i|
			vec[i] = self[i,i]
		}
		return vec
	end

	#=======================================================================
	# type conversion
	#=======================================================================
	
	def to_a
		#: Converts matrix to Array of Arrays
		#| Matrix[[1,2],[3,4]].to_a
		#| => [[1.0, 2.0], [3.0, 4.0]]
		dataclone = []
		@matrix_data_struct.each {|row| dataclone << row.dup}
		return dataclone
	end

	def to_vector
		#: Converts column matrix to Vector. Raise exception if the
		#: matrix is wider than 1.
		raise(DimensionError,"Matrix.to_vector works only with single column matrices") if n != 1
		return Vector.from_array(self.to_a.flatten)
	end

	#=======================================================================
	# iterators
	#=======================================================================
	
	#!
	#=======================================================================
	# Checks
	#=======================================================================
	
	def raise_if_not_equal_size(m2)
		raise(DimensionError,"Matrices must be of equal size") if m != m2.m || n != m2.n
	end
	private :raise_if_not_equal_size

	def raise_if_not_diagonal
		raise(DimensionError,"Matrix must be diagonal") if m != n 
	end
	private :raise_if_not_diagonal
	#!!

	#=======================================================================
	# operators
	#=======================================================================
	
	def +(matrix)
		#: matrix addition
		raise(TypeError, "Argument is not a matrix") unless matrix.kind_of?(Matrix)
		raise_if_not_equal_size(matrix)
		a = []
		m.times{|row|
			a << []
			n.times{|col|
				a[row][col] = @matrix_data_struct[row][col] + matrix[row,col]
			}
		}
		return self.class.new(a)
	end

	def -(matrix)
		raise(TypeError, "Argument is not a matrix") unless matrix.kind_of?(Matrix)
		#: matrix subtraction
		raise_if_not_equal_size(matrix)
		a = []
		m.times{|row|
			a << []
			n.times{|col|
				a[row][col] = @matrix_data_struct[row][col] - matrix[row,col]
			}
		}
		return self.class.new(a)
	end

	def plus!(matrix)
		#: In-place matrix addition, faster than using a = a + b
		#| old_id = a.object_id
		#| a.plus! b
		#| a.object_id == old_id
		#| => true
		raise(TypeError, "Argument is not a matrix") unless matrix.kind_of?(Matrix)
		raise_if_not_equal_size(matrix)
		a = []
		m.times{|row|
			a << []
			n.times{|col|
				@matrix_data_struct[row][col] += matrix[row,col]
			}
		}
		return nil
	end

	def minus!(matrix)
		#: In-place matrix subtraction, faster than using a = a - b<br />
		#: See -->plus! for details
		raise(TypeError, "Argument is not a matrix") unless matrix.kind_of?(Matrix)
		raise_if_not_equal_size(matrix)
		a = []
		m.times{|row|
			a << []
			n.times{|col|
				@matrix_data_struct[row][col] -= matrix[row,col]
			}
		}
		return nil
	end

	def mul!(operand)
		#: In-place multiplication by number
		raise(TypeError,"mul! can only multiply matrix by number") unless operand.kind_of?(Numeric)
		m.times{|row|
			n.times{|col|
				@matrix_data_struct[row][col] *= operand
			}
		}
		return nil
	end

	def matrix_multiply(operand)
		# check
		raise(TypeError, "Argument is not a matrix") unless operand.kind_of?(Matrix)
		raise(DimensionError,"Cant multiply matrices with wrong size") unless n == operand.m
		a = []
		m.times{|row|
			a << []
			operand.n.times{|col|
				x = 0.0
				# = row * operand.col
				n.times {|n| x += self[row,n] * operand[n,col]}
				a[row][col] = x
			}
		}
		return self.class.new(a)
	end

	def matrix_vector_multiply(operand)
		raise(TypeError, "Argument is not a vector") unless operand.kind_of?(Vector)
		raise(DimensionError,"Cant multiply matrix with wector of wrong size") unless n == operand.size
		a = []
		m.times{|i|
			sum = 0.0
			n.times{|j|
				sum += self[i,j]*operand[j]
			}
			a[i] = sum
		}
		return operand.class.from_array(a)
	end

	def comutative_multiply(operand)
		#: Multiplication of matrix by: Float, Fixnum
		if operand.kind_of?(Numeric)
			a = []
			m.times{|row|
				a << []
				n.times{|col|
					a[row][col] = @matrix_data_struct[row][col] * operand
				}
			}
			return self.class.new(a)
		else
			raise(TypeError,"Matrix.*: Can't comutatively multiply matrix by class #{operand.class}")
		end
	end

	def multiply_yield_diagonal(operand)
		# check
		raise(TypeError, "Argument is not a matrix") unless operand.kind_of?(Matrix)
		raise(DimensionError,"Cant multiply matrices with wrong size") unless n == operand.m
		size = [self.m, operand.n].min
		a = []

		size.times{|i|
			x = 0.0
			# = row * operand.col
			n.times {|n| x += self[i,n] * operand[n,i]}
			a[i] = x
		}
		return Vector.new(a)
	end

	#=======================================================================
	# comparison
	#=======================================================================
	
	def == (matrix)
		#: Exact comparison of two matrices, element by element
		return false unless matrix.kind_of?(Matrix)
		raise_if_not_equal_size(matrix)
		m.times{|row|
			n.times{|col|
				return false if @matrix_data_struct[row][col] != matrix[row,col]
			}
		}
		return true
	end

	def =~ (matrix)
		#: Loose comparison of two matrices, element by element. Max. allowed
		#: error is stored in class variable, accessible by methods
		#: -->self.epsilon and -->self.epsilon= .
		#| Matrix.epsilon = 0.001
		#| Matrix[[1,2],[3,4]] == Matrix[[1,2],[3,4.0001]]
		#| => false
		#| Matrix[[1,2],[3,4]] =~ Matrix[[1,2],[3,4.0001]]
		#| => true
		return false unless matrix.kind_of?(Matrix)
		raise_if_not_equal_size(matrix)
		m.times{|row|
			n.times{|col|
				return false if (@matrix_data_struct[row][col] - matrix[row,col]).abs > @@default_epsilon
			}
		}
		return true
	end

	#=======================================================================
	# Trace
	#=======================================================================
	
	def trace
		#: Returns trace, the sum of diagonal elements in the matrix
		raise_if_not_diagonal
		sum = 0.0
		m.times{|i|
			sum += @matrix_data_struct[i][i]
		}
		return sum
	end
	
	#=======================================================================
	# Transpose
	#=======================================================================
	
	def transpose
		#: Returns transposed matrix
		a = []
		n.times{|col|
			a << []
			m.times{|row|
				a[col][row] = @matrix_data_struct[row][col]
			}
		}
		return self.class.new(a)
	end

	def transpose!
		#: In-place transpose
		a = []
		n.times{|col|
			a << []
			m.times{|row|
				a[col][row] = @matrix_data_struct[row][col]
			}
		}
		@matrix_data_struct = a
		return nil
	end

end
end

