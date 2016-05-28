module Algebra

# Exceptions
class  SparseWriteError < RuntimeError
end

class  SparseTypeError < RuntimeError
end

class SparseMatrixCSX

	#=======================================================================
	# Constructor
	#=======================================================================

	def initialize(row_mode, m, n, values, indices, record_i)
		# Row mode
		case row_mode
		when :row
			@row_mode = true
		when :col
			@row_mode = false
		else
			raise(SparseTypeError, "Row mode must be either :row or :col")
		end

		@m = m
		@n = n

		if @row_mode
			# Size checks, CSR
			raise(RuntimeError, "Size of row indices array must be m+1") unless record_i.size == m + 1
			raise(RuntimeError, "Size of values and column indices arrays must be the same") unless values.size == indices.size
		else
			# Size checks, CSC
			raise(RuntimeError, "Size of col indices array must be n+1") unless record_i.size == n + 1
			raise(RuntimeError, "Size of values and row indices arrays must be the same") unless values.size == indices.size
		end

		@data_values = values
		@data_indices = indices
		@data_records = record_i
	end

	#=======================================================================
	# Basics
	#=======================================================================
	
	def deep_copy
		if @row_mode
			return SparseMatrixCSX.new(:row, @m, @n, @data_values.dup, @data_indices.dup, @data_records.dup)
		else
			return SparseMatrixCSX.new(:col, @m, @n, @data_values.dup, @data_indices.dup, @data_records.dup)
		end
	end

	def rewrite(newtype)
		if newtype == type
			# Same type - just copy
			return deep_copy
		else
			# Rewrite matrix to new type
			raise_slow
			values = []
			indices = []
			records = []
			if newtype == :row
				count = 0
				m.times{|i|
					records << count
					row_each_nonzero_with_index(i){|v,j|
						values[count] = v
						indices[count] = j
						count += 1

					}
				}
				records << count
			elsif newtype == :col
				count = 0
				n.times{|j|
					records << count
					col_each_nonzero_with_index(j){|v,i|
						values[count] = v
						indices[count] = i
						count += 1
					}
				}
				records << count
			else
				raise(SparseTypeError, "Row mode must be either :row or :col")
			end
			return self.class.new(newtype, m, n, values, indices, records)
		end
	end

	#=======================================================================
	# Matrix information
	#=======================================================================

	def type
		return :row if @row_mode
		return :col
	end

	def change_type!(newtype)
		oldtype = type
		# Assign new type
		case newtype
		when :row
			@row_mode = true
		when :col
			@row_mode = false
		else
			raise(SparseTypeError, "Row mode must be either :row or :col")
		end
		# Swap size info if type changes
		if newtype != oldtype
			@m, @n = @n, @m
		end
		return nil
	end

	def m
		return @m
	end

	def n
		return @n
	end

	#=======================================================================
	# Conversion to Matrix and SparseMatrix
	#=======================================================================

	def to_matrix
		mat = Matrix.zero(m, n)
		m.times{|i|
			row_each_nonzero_with_index(i){|x,j|
				 mat[i,j] = x
			}
		}
		return mat
	end

	def to_sparse_matrix
		mat = SparseMatrix.new(m, n)
		m.times{|i|
			row_each_nonzero_with_index(i){|x,j|
				 mat[i,j] = x
			}
		}
		return mat
	end

	#=======================================================================
	# Element access
	#=======================================================================
	
	def [](i,j)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		raise(IndexError,"Col index (j) out of bounds") if j >= @n || j < 0

		# Indexes in the internal record structure
		if @row_mode
			ri = i
			rj = j
		else
			ri = j
			rj = i
		end

		# Record start and end
		rs = @data_records[ri];
		re = @data_records[ri+1];
		# Traverse the row
		for c in rs...re do
			return @data_values[c] if @data_indices[c] == rj
		end
		# Element not found -> return zero
		return 0.0
	end

	def []=(i,j,value)
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0
		raise(IndexError,"Col index (j) out of bounds") if j >= @n || j < 0

		# Indexes in the internal record structure
		if @row_mode
			ri = i
			rj = j
		else
			ri = j
			rj = i
		end

		# Row start and end
		rs = @data_records[ri];
		re = @data_records[ri+1];
		# Traverse the row
		for c in rs...re do
			if @data_indices[c] == rj
				@data_values[c] = value
				return value
			end

		end
		# Element not found -> raise error
		raise(SparseWriteError, "SparseMatrixCSX can not write into element which is zero")
	end


	#=======================================================================
	# Row / col access
	#=======================================================================
	
	def row_as_array(i)
		array = Array.new(@n, 0.0)
		row_each_nonzero_with_index(i){|x,j|
			array[j] = x
		}
		return array
	end

	def col_as_array(i)
		array = Array.new(@m, 0.0)
		col_each_nonzero_with_index(i){|x,j|
			array[j] = x
		}
		return array
	end

	#=======================================================================
	# Iterators
	#=======================================================================
	
	def row_each_nonzero_with_index(i)
		# Bounds check
		raise(IndexError,"Row index (i) out of bounds") if i >= @m || i < 0

		if @row_mode
			# In row mode, it is easy
			# Row start and end
			rs = @data_records[i];
			re = @data_records[i+1];
			# Traverse the row
			for c in rs...re do
				yield(@data_values[c], @data_indices[c])
			end
		else
			# Slow!
			raise_slow
			@n.times{|j|
				if (val = self[i,j]) != 0
					yield(val,j)
				end
			}
		end
		return nil
	end

	def col_each_nonzero_with_index(j)
		# Bounds check
		raise(IndexError,"Col index (j) out of bounds") if j >= @n || j < 0

		if @row_mode
			# Slow!
			raise_slow
			@m.times{|i|
				if (val = self[i,j]) != 0
					yield(val,i)
				end
			}
		else
			# In col mode, it is easy
			# Row start and end
			rs = @data_records[j];
			re = @data_records[j+1];
			# Traverse the row
			for c in rs...re do
				yield(@data_values[c], @data_indices[c])
			end
		end
		return nil
	end

	#=======================================================================
	# In-matrix operations
	#=======================================================================
	
	def record_dot_record(i1, i2)
		# Record starts and ends
		rs1 = @data_records[i1];
		re1 = @data_records[i1+1];

		rs2 = @data_records[i2];
		re2 = @data_records[i2+1];

		sum = 0.0

		# Traverse the rows simultaneously, using two pointing indices
		p1 = rs1
		p2 = rs2

		while p1 < re1 && p2 < re2 do
			di1 = @data_indices[p1]
			di2 = @data_indices[p2]
			if di1 == di2
				sum += @data_values[p1] * @data_values[p2]
				p1 += 1
				p2 += 1
			elsif di1 < di2
				p1 += 1
			else # di1 > di2
				p2 += 1
			end
		end

		return sum
	end
	

	#=======================================================================
	# Solve
	#=======================================================================

	#  not implemented in ruby version	
end
end
