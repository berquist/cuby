################################################################################
#
#  class Matrix : Common methods
#
################################################################################

#: Methods used by both ruby and binary version of Matrix class

### Matrix-vector multiplication
### Iterators - fix access

module Algebra
require @@directory+"/matrix_common_universal.rb"
require @@directory+"/matrix_plot.rb"

class Matrix
	include MatrixCommonUniversal
	include MatrixPlot
	# Ithis adds methods: to_s, inspect, to_s_rb

	DEFAULT_EPSILON = 1.0e-8
	@@default_epsilon = DEFAULT_EPSILON

	#=======================================================================
	# Constructors
	#=======================================================================

	def self.[](*array)
		#: Creates matrix from array of arrays (used without brackets)
		#| a = Matrix[[1,2],[3,4]]
		#| =>
		#| |1.0 2.0|
		#| |3.0 4.0|
		return from_array(array)
	end

	def self.rows(array)
		#: Synonym for Matrix.from_array
		return from_array(array)
	end

	def self.columns(array)
		#: Creates matrix by columns from array of arrays
		#| a = Matrix([[1,2],[3,4]])
		#| =>
		#| |1.0 3.0|
		#| |2.0 4.0|
		return from_array(array).transpose
	end

	#> self.diagonal(size,diagonal_element)
	#> self.diagonal(array)
	#> self.diagonal(vector)
	#!
	def self.diagonal(arg1, arg2 = 1.0)
	#!!
		#: Creates diagonal matrix
		#* from size and repeated diagonal element (1.0 if ommitted): Matrix.diagonal(size,number)
		#* using array in diagonal: Matrix.diagonal(array)

		case arg1
		when Fixnum
			return self.diagonal1(arg1,arg2)
		when Array
			return self.diagonal2(arg1)
		when Vector
			return self.diagonal2(arg1.to_a)
		else
			raise "Matrix.diagonal: wrong type of arguments"
		end
	end

	def self.identity(n)
		#: Creates identity matrix of size n x n
		return self.diagonal(n,1.0)
	end

	def self.unit(n)
		#: Creates identity matrix of size n x n
		return self.diagonal(n,1.0)
	end

	### def self.join_column_vectors(vectors)
	### def self.join_row_vectors(vectors)
	
	def self.random(m,n = nil)
		#: Creates a new matrix filled with random values 0.0 <= r < 1.0
		x = self.zero(m,n)
		x.each_index{|i,j|
			x[i,j] = rand
		}
		return x
	end

	#=======================================================================
	# basics
	#=======================================================================

	def clone
		#: Synonym for -->deep_copy
		return deep_copy
	end

	#=======================================================================
	# element access
	#=======================================================================
	
	def paste!(i,j,matrix)
		#: Copies smaller matrix to given position at large one
		matrix.m.times{|ii|
			matrix.n.times{|jj|
				self[ii+i,jj+j] = matrix[ii,jj]
			}
		}
		return nil
	end

	def swap_columns!(a,b) # => nil
		# Swap two columns of a matrix
		return nil if a == b
		tmp = nil
		m.times{|i|
			tmp = self[i,a]
			self[i,a] = self[i,b]
			self[i,b] = tmp
		}
		return nil
	end

	def reorder_columns!(array_of_indices)
		#: Reorder colum in matrices according to array of indices.
		n.times{|j|
			k = array_of_indices[j]
			swap_columns!(j,k)
			array_of_indices.each_index{|ii|
				if array_of_indices[ii] == j
					array_of_indices[ii] = k
				elsif array_of_indices[ii] == k
					array_of_indices[ii] = j
				end
			}
		}
		return nil
	end

	#=======================================================================
	# iterators
	#=======================================================================
	
	#> each {|e| ... }
	#: Iterate over all elements
	#| Matrix[[1,2],[3,4]].each {|e| puts e}
	#|
	#| 1.0
	#| 3.0
	#| 2.0
	#| 4.0
	#!
	def each
	#!!
		m.times{|row|
			n.times{|col|
				yield(self[row,col])
			}
		}
	end

	#> each_index {|i,j| ... }
	#: Iterate over element indexes, i = index in row, j = index of column
	#| Matrix[[1,2],[3,4]].each_index {|i,j| puts "#{i},#{j}"}
	#|
	#| 0,0
	#| 1,0
	#| 0,1
	#| 1,1
	#!
	def each_index
	#!!
		n.times{|col|
			m.times{|row|
				yield(row,col)
			}
		}
	end

	#> each_with_index {|e,i,j| ... }
	#: Iterate over all elements, yielding both the elemnt and its index<br />
	#| Matrix[[1,2],[3,4]].each_with_index {|e,i,j| puts "#{i},#{j} -> #{e}"}
	#|
	#| 0,0 -> 1.0
	#| 1,0 -> 3.0
	#| 0,1 -> 2.0
	#| 1,1 -> 4.0
	#!
	def each_with_index
	#!!
		n.times{|col|
			m.times{|row|
				yield(self[row,col],row,col)
			}
		}
	end

	def row_each_with_index(row)
	#!!
		n.times{|col|
			yield(self[row,col],col)
		}
	end

	def each_diagonal
		[m,n].min.times {|i|
			yield self[i,i]
		}
	end

	### map
	### map_with_index
	### map!
	### map_with_index!

	#!

	#=======================================================================
	# Norms etc.
	#=======================================================================
	
	def min
		#: Return the smallest number in the matrix
		minimum = self[0,0]
		each{|x| minimum = x if x < minimum}
		return minimum
	end
	
	def max
		#: Return the largest number in the matrix
		maximum = self[0,0]
		each{|x| maximum = x if x > maximum}
		return maximum
	end
	
	def min_with_index
		#: Return the smallest number in the matrix and indices of it 
		minimum = self[0,0]
		i = 0
		j = 0
		each_with_index{|x, row, col| 
			if x < minimum
				minimum = x 
				i = row
				j = col
			end
		}
		return [minimum, i, j]
	end

	def max_with_index
		#: Return the largest number in the matrix and indices of it 
		maximum = self[0,0]
		i = 0
		j = 0
		each_with_index{|x, row, col| 
			if x > maximum
				maximum = x 
				i = row
				j = col
			end
		}
		return [maximum, i, j]
	end

	def count_nonzero(threshold = nil)
		count = 0
		each{|element|
			if threshold
				count += 1 if element.abs > threshold
			else
				count += 1 if element != 0
			end
		}
		return count
	end

	#=======================================================================
	# unary operators
	#=======================================================================
	
	def +@
		#: Unary plus. Does nothing (returns self)
		return self
	end

	def -@
		#: Unary minus. Equal to * -1.0
		return self * -1.0
	end
	
	#=======================================================================
	# operators
	#=======================================================================
	

	def div!(operand)
		#: In-place division
		return self.mul!(1.0/operand)
	end

	#=======================================================================
	# operators: matrix - something
	#=======================================================================
	
	def *(operand)
		#: Multiplication of matrix by: Float, Fixnum, Vector and Matrix
		if operand.kind_of?(Numeric)
			return self.comutative_multiply(operand)
		elsif operand.kind_of?(Matrix)
			return self.matrix_multiply(operand)
		elsif operand.kind_of?(Vector)
			return self.matrix_vector_multiply(operand)
		else
			raise(TypeError,"Matrix.*: Can't multiply matrix by class #{operand.class}")
		end
	end

	def /(operand)
		#: Division of matrix by scalar
		return self.comutative_multiply(1.0/operand)
	end

	def vt_self_v(v1, v2)
		#: Shortcut (V1.to_matrix.transpose * self * v2)
		#!# Can be implemented more efficiently
		raise(ArgumentError, "Argument v1 must be Vector") unless v1.kind_of?(Vector)
		raise(ArgumentError, "Argument v2 must be Vector") unless v2.kind_of?(Vector)
		return (v1.to_matrix.transpose * self * v2)[0]
	end

	#=======================================================================
	# comparison
	#=======================================================================
	
	def self.epsilon
		#: Returns threshold for loose comparison using -->=~ operator.
		#: The value is shared between Matrix and Vector classes
		return @@default_epsilon
	end

	def self.epsilon=(value)
		#: Sets threshold for loose comparison using -->=~ operator.
		#: The value is shared between Matrix and Vector classes
		@@default_epsilon = value
	end

	#=======================================================================
	# Matrix operations
	#=======================================================================
	
	def inverse2
		a = self[0,0]
		b = self[0,1]
		c = self[1,0]
		d = self[1,1]

		return Matrix[[d, -b],[-c, a]] * (1.0 / (a*d - b*c))
	end

	def inverse3
		a = self[0,0]
		b = self[0,1]
		c = self[0,2]
		d = self[1,0]
		e = self[1,1]
		f = self[1,2]
		g = self[2,0]
		h = self[2,1]
		k = self[2,2]

		det = a*(e*k - f*h) + b*(f*g - k*d) + c*(d*h - e*g)
		if det == 0.0
			raise "Matrix is not invertible"
		end

		m = Matrix[
			[e*k-f*h,c*h-b*k,b*f-c*e],
			[f*g-d*k,a*k-c*g,c*d-a*f],
			[d*h-e*g,b*g-a*h,a*e-b*d]
		]

		return m * (1.0/det)
	end

	#=======================================================================
	# Matrix operations - shortcuts
	#=======================================================================
	
	def det
		#: Synonym for -->determinant
		return determinant
	end

	def inv
		#: Shortcut for -->inverse
		inverse
	end
	
	def inv!
		#: Shortcut for -->inverse!
		inverse!
	end

	def pseudoinverse(threshold = 1.0e-10)
		#: Moore-penrose pseudoinverse using SVD decomposition.
		#: Threshold parameter is used to evaluate zero/nonzero
		#: elements.
		u, s, v = self.svd
		[s.m, s.n].min.times {|i|
			s[i,i] = 1.0 / s[i,i] if s[i,i].abs > threshold
		}
		return v.transpose * s.transpose * u.transpose
	end

	#=======================================================================
	# Eigensystem sorter
	#=======================================================================
	
	def Matrix.sort_eigensystem!(vec, real, imag, sort = :descending_real) # => nil
		#: Sort eigenvalues and eigenvectors (original matrices are modified). Possible sorting:
		#* :ascending_real
		#* :descending_real
		#* :ascending_real_abs
		#* :descending_real_abs

		# Build array of elements to be sorted: [real, imag, vector_column_index]
		array = []
		vec.n.times{|i|
			array[i] = [real[i,0], imag[i,0], i]
		}

		# Sort
		case sort
		when :ascending_real
			array.sort!{|x,y| x[0] <=> y[0] }
		when :descending_real
			array.sort!{|x,y| y[0] <=> x[0] }
		when :ascending_real_abs
			array.sort!{|x,y| x[0].abs <=> y[0].abs }
		when :descending_real_abs
			array.sort!{|x,y| y[0].abs <=> x[0].abs }
		else
		end

		# Rebuild matrices
		reorder_prescription = []
		vec.n.times{|i|
			real[i,0] = array[i][0]
			imag[i,0] = array[i][1]
			reorder_prescription[i] = array[i][2]
		}
		vec.reorder_columns!(reorder_prescription)

		return reorder_prescription
	end
	
	#=======================================================================
	# File read / write
	#=======================================================================
	
	def to_file(file)
		#: Write matrix to text file. File can be specified
		#: as IO object or filename.
		close = false
		if file.class == String
			close = true
			file = File.new(file,"w+")
		end

		m.times {|i|
			n.times{|j|
				file.printf("%16.8e",self[i,j])	
			}
			file.puts
		}

		file.close if close
		return nil
	end

	def Matrix.from_file(file)
		#: Read matrix from text file. File can be specified
		#: as IO object or filename.
		close = false
		if file.class == String
			close = true
			file = File.new(file,"r")
		end

		a = []
		file.each_line {|line|
			if line !~ /^\s*$/
				a << line.split.map{|x| x.to_f}
			end
		}

		new_matrix = Matrix.from_array(a)

		file.close if close
		return new_matrix
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
					arr = self.to_a
					arr.each{|a| def a.to_yaml_style; :inline; end}
					map.add('elements', arr)
				end
			end
		end

		def self.yaml_new(klass, tag, val)
			array = val['elements']
			return Matrix.from_array(array)
		end

	elsif yaml_engine == 'psych'
		# Psych yaml engine
	
		def encode_with(coder)
			coder.tag = "!cuby.molecular.cz,2009/Algebra::Matrix"
			arr = self.to_a
			coder.style = Psych::Nodes::Mapping::BLOCK # BLOCK or FLOW
			#!# Would be nice to have rows styled as FLOW, but this will be more complicated
			coder['elements'] = arr
		end

		def init_with(coder)
			puts "init_with"
		end

		YAML::add_domain_type('cuby.molecular.cz,2009', 'Algebra::Matrix') { |type, value|
			Matrix.from_array(value['elements'])
		}

	else
		raise "YAML engine can not be determined"
	end
end
end

