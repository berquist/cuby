module Algebra
require @@directory+"/matrix_common_universal.rb"

class  SparseSlowError < RuntimeError
end

class SparseMatrixCSX

	include MatrixCommonUniversal
	# This adds methods: to_s, inspect, to_s_rb
	
	#=======================================================================
	# Constructors
	#=======================================================================

	def self.from_sparse_matrix(row_mode, mat)
		if row_mode == :row
			val = []
			col_i = []
			row_i = []

			count = 0
			mat.m.times{|i|
				row_i << count
				mat.row_nonzero_indices_array(i).sort.each{|j|
					val[count] = mat[i,j]
					col_i[count] = j
					count += 1
				}
			}
			row_i << count
			return self.new(row_mode, mat.m, mat.n, val, col_i, row_i)
		elsif row_mode == :col
			val = []
			data_i = []
			col_i = []
			count = 0

			mat.n.times{|j|
				col_i << count
				mat.m.times{|i|
					a = mat[i,j]
					if a != 0.0
						val[count] = a
						data_i[count] = i
						count += 1
					end
				}
			}
			col_i << count
			return self.new(row_mode, mat.m, mat.n, val, data_i, col_i)
		else
			raise(ArgumentError, "Row mode must be either :row or :col, is #{row_mode}")
		end
	end

	def self.from_matrix(row_mode, mat)
		if row_mode == :row
			val = []
			col_i = []
			row_i = []

			count = 0
			mat.m.times{|i|
				row_i << count
				mat.row_each_with_index(i){|a,j|
					if a != 0.0
						val[count] = a
						col_i[count] = j
						count += 1
					end
				}
			}
			row_i << count
			return self.new(row_mode, mat.m, mat.n, val, col_i, row_i)
		elsif row_mode == :col
			val = []
			data_i = []
			col_i = []
			count = 0

			mat.n.times{|j|
				col_i << count
				mat.m.times{|i|
					a = mat[i,j]
					if a != 0.0
						val[count] = a
						data_i[count] = i
						count += 1
					end
				}
			}
			col_i << count
			return self.new(row_mode, mat.m, mat.n, val, data_i, col_i)
		else
			raise(ArgumentError, "Row mode must be either :row or :col, is #{row_mode}")
		end
	end

	#=======================================================================
	# Slow operations
	#=======================================================================

	@@allow_slow_operations = false

	def raise_slow
		unless @@allow_slow_operations
			raise(SparseSlowError, "SparseMatrixCSX 'slow' operations are disabled")
		end
	end

	def self.allow_slow_operations=(value)
		@@allow_slow_operations = value
	end

	def allow_slow_operations
		return @@allow_slow_operations
	end
	
	#=======================================================================
	# Basics
	#=======================================================================

	def clone
		#: Synonym for -->deep_copy
		return deep_copy
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
	# Algebra
	#=======================================================================
	
	def invert_diag_3x3_blocks!
		#: Inverts all 3x3 blocks along the diagonal. If the rest of teh matrix is zero, this
		#: yield the inverse of teh whole matrix.
		unless self.m % 3 == 0 && self.n == self.m
			raise "Method invert_diag_3x3_blocks! is applicable only to square matrices with size divisible by 3"
		end

		(self.m/3).times{|i|
			bi = i*3
			a = self[bi+0, bi+0]
			b = self[bi+0, bi+1]
			c = self[bi+0, bi+2]
			d = self[bi+1, bi+0]
			e = self[bi+1, bi+1]
			f = self[bi+1, bi+2]
			g = self[bi+2, bi+0]
			h = self[bi+2, bi+1]
			k = self[bi+2, bi+2]

			det = a*(e*k - f*h) + b*(f*g - k*d) + c*(d*h - e*g)
			if det == 0.0
				raise "Matrix is not invertible, determinant of block #{bi+1} is zero"
			end

			self[bi+0, bi+0] = (e*k-f*h)/det
			self[bi+0, bi+1] = (c*h-b*k)/det
			self[bi+0, bi+2] = (b*f-c*e)/det
			self[bi+1, bi+0] = (f*g-d*k)/det
			self[bi+1, bi+1] = (a*k-c*g)/det
			self[bi+1, bi+2] = (c*d-a*f)/det
			self[bi+2, bi+0] = (d*h-e*g)/det
			self[bi+2, bi+1] = (b*g-a*h)/det
			self[bi+2, bi+2] = (a*e-b*d)/det
		}

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

	#--------------------------------------------------
	# def comutative_multiply(operand)
	# 	#: Multiplication of matrix by: Float, Fixnum
	# 	if operand.kind_of?(Numeric)
	# 		result = SparseMatrix.new(m,n)
	# 		each_nonzero_with_index{|x,i,j| result[i,j] = x * operand}
	# 		return result
	# 	else
	# 		raise(TypeError,"Matrix.*: Can't comutatively multiply matrix by class #{operand.class}")
	# 	end
	# end
	#-------------------------------------------------- 

	#=======================================================================
	# Operators: matrix - something
	#=======================================================================
	
	def *(operand)
		#: Multiplication of matrix by: Float, Fixnum, Vector and Matrix
		if operand.kind_of?(Vector)
			if self.respond_to?(:matrix_vector_multiply_c) && self.type == :row
				return self.matrix_vector_multiply_c(operand)
			else
				return self.matrix_vector_multiply(operand)
			end
		#--------------------------------------------------
		# elsif operand.kind_of?(Numeric)
		# 	return self.comutative_multiply(operand)
		#-------------------------------------------------- 
		#--------------------------------------------------
		# elsif operand.kind_of?(Matrix)
		# 	return self.matrix_multiply(operand)
		#-------------------------------------------------- 
		else
			raise(TypeError,"#{self.class}.*: Can't multiply matrix by class #{operand.class}")
		end
	end

	#=======================================================================
	# In-matrix operations
	#=======================================================================
	
	def col_dot_col(j1, j2)
		raise(IndexError,"Col index (j1) out of bounds") if j1 >= m || j1 < 0
		raise(IndexError,"Col index (j2) out of bounds") if j2 >= m || j2 < 0
		if type == :row
			raise_slow 
			s = 0.0
			col_each_nonzero_with_index(j1){|x,k|
				s += x * self[k,j2]
			}
			return s
		else
			return record_dot_record(j1, j2)
		end
	end
	
	def row_dot_row(i1, i2)
		raise(IndexError,"Row index (i1) out of bounds") if i1 >= m || i1 < 0
		raise(IndexError,"Row index (i2) out of bounds") if i2 >= m || i2 < 0
		if type == :col
			raise_slow 
			s = 0.0
			row_each_nonzero_with_index(i1){|x,k|
				s += x * self[i2,k]
			}
			return s
		else
			return record_dot_record(i1, i2)
		end
	end
	
end
end
