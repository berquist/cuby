# Simple implementation of sparse matrices


module Algebra
require @@directory+"/matrix_common_universal.rb"
class SparseMatrix
	include MatrixCommonUniversal
	# This adds methods: to_s, inspect, to_s_rb

	attr_reader :m, :n

	#=======================================================================
	# Constructors
	#=======================================================================
	
	def initialize(m, n = nil)
		# Creates an empty matrix
		raise(TypeError,"Matrix dimension m must be larger than 0") if m <= 0
		n = m if n.nil? # Only one parameter -> square matrix
		raise(TypeError,"Matrix dimension n must be larger than 0") if n <= 0
		@m = m
		@n = n
		@rows = {}
	end

	def self.from_matrix(matrix)
		raise(TypeError,"Argument of SparseMatrix.from_matrix must be Matrix") unless matrix.kind_of?(Matrix)
		result = SparseMatrix.new(matrix.m,matrix.n)
		matrix.each_with_index{|x,i,j|
			result[i,j] = x unless x == 0.0
		}
		return result
	end

	#=======================================================================
	# Element access
	#=======================================================================

	def []=(i,j,value)
		# Check i,j
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		raise(IndexError,"Col index (j) out of bounds") if j >= @n || j < 0
		raise(TypeError,"Only numbers can be stored in matrix") unless value.kind_of?(Numeric)
		value = value.to_f unless value.kind_of?(Float)
		# Erase element instead of setting it to zero
		if value == 0.0
			return 0.0 unless @rows[i] # Return if zero was to be written into new row
			return 0.0 unless @rows[i].has_key?(j)
			@rows[i].delete(j) # Remove element j
			@rows.delete(i) if @rows[i].empty? # Remove row when it stores no values
			return 0.0
		end
		# Create new row if needed
		# It is a hash that returns 0.0 if the key does not exist
		@rows[i] = Hash.new(0.0) unless @rows[i]
		# Save value
		@rows[i][j] = value
	end

	def [](i,j)
		# Check i,j
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		raise(IndexError,"Col index (j) out of bounds") if j >= @n || j < 0
		# Zero if row does not exist
		return 0.0 unless @rows[i]
		return @rows[i][j]
	end

	def nonzero?(i,j)
		return false unless @rows[i]
		return false unless @rows[i].has_key?(j)
		return true
	end

	def nonzero_count
		count = 0
		@rows.each_value{|row| count += row.size}
		return count
	end

	#=======================================================================
	# Row access
	#=======================================================================

	def row_nonzero_count(i)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		return 0 unless @rows[i]
		return @rows[i].size
	end

	def row_as_array(i)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		a = Array.new(n,0.0)
		row_each_nonzero_with_index(i){|x,j| a[j] = x}
		return a
	end

	def row_nonzero_as_array(i)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		return [] unless @rows[i]
		return @rows[i].values
	end

	def row_nonzero_indices_array(i)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		return [] unless @rows[i]
		return @rows[i].keys
	end

	#=======================================================================
	# Row iterators
	#=======================================================================

	def row_as_vector(i)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		a = Vector.zero(n)
		row_each_nonzero_with_index(i){|x,j| a[j] = x}
		return a
	end

	def row_each_nonzero_with_index(i)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		return nil unless @rows[i]
		@rows[i].each_pair{|j, value|
			yield(value,j)
		}
		return nil
	end

	def row_each_nonzero_index(i)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		return nil unless @rows[i]
		@rows[i].each_key{|j|
			yield(j)
		}
		return nil
	end

	#=======================================================================
	# Iterators
	#=======================================================================

	def each_nonzero_with_index
		@rows.each_pair{|i,row|
			row.each_pair{|j,value|
				yield(value,i,j)
			}
		}
		return nil
	end

	def each_nonzero_index
		@rows.each_pair{|i,row|
			row.each_key{|j|
				yield(i,j)
			}
		}
		return nil
	end

	#=======================================================================
	# Operators
	#=======================================================================
	
	def matrix_vector_multiply(operand)
		raise(TypeError, "Argument is not a vector") if operand.class != Vector
		raise(DimensionError,"Cant multiply matrix with vector of wrong size") unless n == operand.size
		a = []
		m.times{|i|
			sum = 0.0
			row_each_nonzero_with_index(i){|x,j|
				sum += x * operand[j]
			}
			a[i] = sum
		}
		return Vector.from_array(a)
	end

	def comutative_multiply(operand)
		#: Multiplication of matrix by: Float, Fixnum
		if operand.kind_of?(Numeric)
			result = SparseMatrix.new(m,n)
			each_nonzero_with_index{|x,i,j| result[i,j] = x * operand}
			return result
		else
			raise(TypeError,"Matrix.*: Can't comutatively multiply matrix by class #{operand.class}")
		end
	end

	#=======================================================================
	# Operators: matrix - something
	#=======================================================================
	
	def *(operand)
		#: Multiplication of matrix by: Float, Fixnum, Vector and Matrix
		if operand.kind_of?(Numeric)
			return self.comutative_multiply(operand)
		#--------------------------------------------------
		# elsif operand.kind_of?(Matrix)
		# 	return self.matrix_multiply(operand)
		#-------------------------------------------------- 
		elsif operand.kind_of?(Vector)
			return self.matrix_vector_multiply(operand)
		else
			raise(TypeError,"SparseMatrix.*: Can't multiply matrix by class #{operand.class}")
		end
	end

	#=======================================================================
	# Type conversion
	#=======================================================================

	def to_matrix
		m = Matrix.zero(@m,@n)
		each_nonzero_with_index{|x,i,j| m[i,j] = x}
		return m
	end

	def to_sparse_matrix_csx(row_mode)
		return SparseMatrixCSX.from_sparse_matrix(row_mode, self)
	end

	#=======================================================================
	# Comparison: not implemented
	#=======================================================================
	
	def ==(x)
		raise(NoMethodError, "Sparse matrix comparison not implemented")
	end

	def =~(x)
		raise(NoMethodError, "Sparse matrix comparison not implemented")
	end

	#=======================================================================
	# Matrix Market Exchange Format
	#=======================================================================
	
	def write_mmexch_sparse(file = $stdout)
		file.puts "%%MatrixMarket matrix coordinate real general"
		# Header: m, n, c
		file.printf("%6d%6d%8d", m, n, nonzero_count)	
		# Data, indexes i and j start from 1
		each_nonzero_with_index{|x,i,j|
			file.printf("\n%6d%6d%18.6e", i+1, j+1, x)
		}
		return nil
	end

end
end
