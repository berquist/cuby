if ARGV[0] == 'ruby'
	$algebra_force_ruby = true
	ARGV.delete_at(0)
end

require "test/unit"
require "fileutils"
require "../algebra"
include Algebra

puts "--------------------------------------------------------------------------------"
puts "Using version: #{Algebra::version}"
puts "--------------------------------------------------------------------------------"

class TestSparseMatrix < Test::Unit::TestCase

	def test_init
		# Constructors, should return instance of SparseMatrix
		assert_instance_of(SparseMatrix, SparseMatrix.new(5))
		assert_instance_of(SparseMatrix, SparseMatrix.new(5,6))
		assert_instance_of(SparseMatrix, SparseMatrix.from_matrix(Matrix.random(7,8)))

		# Check size of the matrix created
		m = SparseMatrix.new(5)
		assert_equal(5,m.m)
		assert_equal(5,m.n)
		m = SparseMatrix.new(10,20)
		assert_equal(10,m.m)
		assert_equal(20,m.n)
		m = SparseMatrix.from_matrix(Matrix[[1,2],[3,4],[5,6]])
		assert_equal(3,m.m)
		assert_equal(2,m.n)

		# Check contents 
		m = SparseMatrix.from_matrix(Matrix[[1,2],[3,4],[5,6]])
		assert_equal(6, m[2,1])
		m = SparseMatrix.new(10,20)
		assert_equal(0, m[2,1])

		# Check error handling in SparseMatrix.from_matrix
		assert_raise(TypeError) { SparseMatrix.from_matrix([])}

		# Zero and negative-sized matrices
		assert_raise(TypeError) { SparseMatrix.new(0) }
		assert_raise(TypeError) { SparseMatrix.new(-1) }
		assert_raise(TypeError) { SparseMatrix.new(1,0) }
		assert_raise(TypeError) { SparseMatrix.new(0,1) }
		assert_raise(TypeError) { SparseMatrix.new(-1,1) }
		assert_raise(TypeError) { SparseMatrix.new(1,-1) }
	end

	def test_element_access
		# Read/write element
		m = SparseMatrix.from_matrix(Matrix[[1,2],[3,4]])
		assert_equal(1,m[0,0])
		assert_equal(2,m[0,1])
		assert_equal(3,m[1,0])
		assert_equal(4,m[1,1])
		m[1,0] = 15
		assert_equal(15,m[1,0])

		# Exceptions on read/write
		assert_raise(IndexError) {m[0,3]}
		assert_raise(IndexError) {m[3,0]}
		assert_raise(IndexError) {m[0,-3]}
		assert_raise(IndexError) {m[-3,0]}
		assert_raise(IndexError) {m[0,3] = 15}
		assert_raise(IndexError) {m[3,0] = 12}
		assert_raise(IndexError) {m[0,-3] = 15}
		assert_raise(IndexError) {m[-3,0] = 12}
		assert_raise(TypeError)  {m[0,0] = "foo"}

		# nonzero?, nonzero_count
		m = SparseMatrix.new(5,5)
		m[1,2] = 0.05;
		m[2,2] = 1.0;
		m[2,4] = -0.05;
		m[3,2] = -1.0;
		m[4,3] = 0.0;

		assert_equal(4,m.nonzero_count)
		assert_equal(true, m.nonzero?(1,2))
		assert_equal(true, m.nonzero?(2,4))
		assert_equal(false, m.nonzero?(4,3))
		assert_equal(false, m.nonzero?(0,3))
	end

	def test_row_access
		m = SparseMatrix.new(5,5)
		m[1,2] = 0.05;
		m[2,2] = 1.0;
		m[2,4] = -0.05;
		m[3,2] = -1.0;
		m[4,3] = 0.0;

		# row_nonzero_count
		assert_equal(0, m.row_nonzero_count(0))
		assert_equal(1, m.row_nonzero_count(1))
		assert_equal(2, m.row_nonzero_count(2))
		assert_equal(1, m.row_nonzero_count(3))
		assert_equal(0, m.row_nonzero_count(4))

		# row_as_array
		assert_equal([0.0, 0.0, 0.0, 0.0, 0.0],m.row_as_array(0));
		assert_equal([0.0, 0.0, 1.0, 0.0, -0.05],m.row_as_array(2));

		# row_nonzero_as_array, row_nonzero_indices_array
		assert_equal([], m.row_nonzero_as_array(4))
		assert_equal([], m.row_nonzero_indices_array(4))

		assert_equal([1.0, -0.05], m.row_nonzero_as_array(2))
		assert_equal([2, 4], m.row_nonzero_indices_array(2))
	end

	def test_row_iterators
		m = SparseMatrix.new(5,5)
		m[1,2] = 0.05;
		m[2,2] = 1.0;
		m[2,4] = -0.05;
		m[3,2] = -1.0;
		m[4,3] = 0.0;

		#row_each_nonzero_with_index
		av = []
		aj = []
		m.row_each_nonzero_with_index(0){|e,j|
			av << e
			aj << j
		}
		assert_equal([], av)
		assert_equal([], aj)

		av = []
		aj = []
		m.row_each_nonzero_with_index(2){|e,j|
			av << e
			aj << j
		}
		assert_equal([1.0, -0.05], av)
		assert_equal([2, 4], aj)

		# change order of writes - elements should be sorted
		m = SparseMatrix.new(5,5)
		m[1,2] = 0.05;
		m[2,4] = -0.05;
		m[2,2] = 1.0;
		m[3,2] = -1.0;
		m[4,0] = 1.0;
		m[4,1] = 2.0;
		m[4,3] = 0.0;

		av = []
		aj = []
		m.row_each_nonzero_with_index(2){|e,j|
			av << e
			aj << j
		}
		#!# assert_equal([1.0, -0.05], av)
		#!# assert_equal([2, 4], aj)

		#row_each_nonzero_index
		aj = []
		m.row_each_nonzero_index(4){|j|
			aj << j
		}
		assert_equal([0,1], aj)
	end

	def test_matrix_something_operators
		# SparseMatrix * Vector
		a = Matrix[[1,2,0],[4,0,0],[0,8,9]]
		m = SparseMatrix.from_matrix(a)
		v = Vector[10,11,12]
		result = Vector[32.0, 40.0, 196.0]

		assert_instance_of(Vector, m.matrix_vector_multiply(v))
		assert_instance_of(Vector, m*v)

		assert_equal(result, m.matrix_vector_multiply(v))
		assert_equal(result, m*v)

		m2 = SparseMatrix.from_matrix(Matrix[[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
		assert_equal(Vector[68.00, 167.00, 266.00, 365.00], m2*v)
		assert_raise(DimensionError) {SparseMatrix.from_matrix(m2.to_matrix.transpose) * v}

		m3 = SparseMatrix.from_matrix(Matrix[[1,2,3],[4,5,6]])
		assert_equal(Vector[68.00, 167.00], m3*v)
		assert_raise(DimensionError) {SparseMatrix.from_matrix(m3.to_matrix.transpose) * v}

		# SparseMatrix * Vector returning single number
		m = SparseMatrix.from_matrix(Matrix[[1,2,3]])
		assert_instance_of(Vector, m * m.to_matrix.transpose.to_vector)
		assert_equal(Vector[14.0], m * m.to_matrix.transpose.to_vector)

		# SparseMatrix * Numeric
		m = SparseMatrix.from_matrix(Matrix[[1,2,0],[4,0,6]])
		r = SparseMatrix.from_matrix(Matrix[[2,4,0],[8,0,12]])

		assert_instance_of(SparseMatrix, m * 2.0)
		# .to_matrix is needed, comparison of sparse matrices not implemented
		assert_equal((m * 2.0).to_matrix, r.to_matrix)
	end

	def test_numeric_extensions
		# Comutative multiply: Numeric * SparseMatrix
		m = SparseMatrix.from_matrix(Matrix[[1,2,0],[4,0,6]])
		r = SparseMatrix.from_matrix(Matrix[[2,4,0],[8,0,12]])

		# Float
		assert_instance_of(SparseMatrix, 2.0 * m)
		# .to_matrix is needed, comparison of sparse matrices not implemented
		assert_equal((2.0 * m).to_matrix, r.to_matrix)

		# Fixnum
		assert_instance_of(SparseMatrix, 2 * m)
		# .to_matrix is needed, comparison of sparse matrices not implemented
		assert_equal((2 * m).to_matrix, r.to_matrix)
	end

	def test_type_conversion
		m = SparseMatrix.new(2)
		assert_instance_of(String, m.to_s)
		assert_instance_of(String, m.inspect)
		assert_instance_of(String, m.to_s_rb)

		assert_instance_of(Matrix, m.to_matrix)
		assert_instance_of(SparseMatrixCSR, m.to_sparse_matrix_csr)
	end

	def test_comparison
		# Should fail, not implemented
		m = SparseMatrix.new(2)
		m2 = SparseMatrix.new(2)

		assert_raise(NoMethodError) {m == m2}
		assert_raise(NoMethodError) {m =~ m2}
	end

	def test_read_write
		# Read not implemented, no way to test
	end

end
