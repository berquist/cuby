if ARGV[0] == 'col'
	$mode = :col
	ARGV.delete_at(0)
elsif ARGV[0] == 'row'
	$mode = :row
	ARGV.delete_at(0)
else
	raise "Mode must be set in the argument, either 'row' or 'col'"
end

if ARGV[0] == 'ruby'
	$algebra_force_ruby = true
	ARGV.delete_at(0)
end

require "test/unit"
require "fileutils"
require "../algebra"
include Algebra

puts "--------------------------------------------------------------------------------"
puts "Using version: #{Algebra::version}, mode: #{$mode}"
puts "--------------------------------------------------------------------------------"

SparseMatrixCSX.allow_slow_operations = true

class TestSparseMatrix < Test::Unit::TestCase


	def TestSparseMatrix.testmatrix
		m = SparseMatrix.new(5,5)
		m[1,2] = 0.05
		m[2,2] = 1.0
		m[2,4] = -0.05
		m[3,2] = -1.0
		m[4,3] = 0.0
		return m.to_sparse_matrix_csx($mode)
	end


	def test_init
		# Constructors, should return instance of SparseMatrix
		assert_instance_of(SparseMatrixCSX, SparseMatrixCSX.from_matrix($mode, Matrix.random(7,8)))
		assert_instance_of(SparseMatrixCSX, SparseMatrixCSX.from_sparse_matrix($mode, SparseMatrix.new(5)))

		# Check size of the matrix created
		m = SparseMatrixCSX.from_sparse_matrix($mode, SparseMatrix.new(2,10))
		assert_equal(m.m, 2)
		assert_equal(m.n, 10)
	end

	def test_basic
		# Dimensions
		m = SparseMatrixCSX.from_matrix($mode, Matrix[[1,0],[3,4],[0,0]])
		assert_equal(3,m.m)
		assert_equal(2,m.n)
		assert_equal([3,2],m.dimensions)
		assert_equal([3,2],m.size)
		assert_equal("3 x 2",m.size_str)

		# Clone / deep_copy
		m = SparseMatrixCSX.from_matrix($mode, Matrix[[1,0],[0,2]])
		m2 = m.clone
		assert_equal(m.to_matrix,m2.to_matrix)
		assert_not_equal(m.object_id,m2.object_id)
		m2[1,1] = 0
		assert_not_equal(m.to_matrix,m2.to_matrix)
	end

	def test_rewrite
		m = TestSparseMatrix.testmatrix

		# Rewrite to the same type == deep copy
		m2 = m.rewrite($mode)
		assert_equal(m.to_matrix,m2.to_matrix)
		assert_not_equal(m.object_id,m2.object_id)
		assert_equal(m.type,m2.type)
		m2[1,2] = 0
		assert_not_equal(m.to_matrix,m2.to_matrix)

		# Get the other mode
		if $mode == :col
			newmode = :row
		else
			newmode = :col
		end

		# Rewrite changing type
		m2 = m.rewrite(newmode)
		assert_not_equal(m.object_id,m2.object_id)
		assert_not_equal(m.type,m2.type)
		assert_equal(m.to_matrix,m2.to_matrix)

		# Type change without modification
		m2 = m.deep_copy
		m2.change_type!(newmode)
		assert_equal(m.to_matrix,m2.to_matrix.transpose)
	end

	def test_element_access
		# Read/write element
		m = SparseMatrix.from_matrix(Matrix[[1,2,0],[3,4,0]]).to_sparse_matrix_csx($mode)
		assert_equal(1,m[0,0])
		assert_equal(2,m[0,1])
		assert_equal(0,m[0,2])
		assert_equal(3,m[1,0])
		assert_equal(4,m[1,1])
		assert_equal(0,m[1,2])

		# Indexing exceptions
		assert_raise(IndexError) {m[0,3]}
		assert_raise(IndexError) {m[3,0]}
		assert_raise(IndexError) {m[0,-3]}
		assert_raise(IndexError) {m[-3,0]}

		# Write
		assert_raise(SparseWriteError) {m[0,2] = 1.0}
		m[0,1] = 15.5
		assert_equal(15.5, m[0,1])
	end

	def test_row_access
		m = TestSparseMatrix.testmatrix

		# row_as_array
		assert_equal([0.0, 0.0, 0.0, 0.0, 0.0],m.row_as_array(0))
		assert_equal([0.0, 0.0, 1.0, 0.0, -0.05],m.row_as_array(2))

		# col as array
		assert_equal([0.0, 0.05, 1.0, -1.0, 0.0], m.col_as_array(2))
		assert_equal([0.0, 0.0, 0.0, 0.0, 0.0], m.col_as_array(3))

	end

	def test_iterators
		m = TestSparseMatrix.testmatrix

		# row - empty
		index = []
		value = []
		m.row_each_nonzero_with_index(0){|v,i| index << i; value << v}
		assert_equal([], index)
		assert_equal([], value)

		# row - with elements
		index = []
		value = []
		m.row_each_nonzero_with_index(2){|v,i| index << i; value << v}
		assert_equal([2, 4], index)
		assert_equal([1.0, -0.05], value)

		# col - empty
		index = []
		value = []
		m.col_each_nonzero_with_index(3){|v,i| index << i; value << v}
		assert_equal([], index)
		assert_equal([], value)

		# col - with elements
		index = []
		value = []
		m.col_each_nonzero_with_index(2){|v,i| index << i; value << v}
		assert_equal([1,2,3], index)
		assert_equal([0.05, 1.0, -1.0], value)
	end

	def test_type_conversion
		m = SparseMatrix.new(2)
		assert_instance_of(SparseMatrixCSX, m.to_sparse_matrix_csx($mode))

		m = SparseMatrixCSX.from_sparse_matrix($mode, m)
		assert_instance_of(String, m.to_s)
		assert_instance_of(String, m.inspect)
		assert_instance_of(String, m.to_s_rb)

		assert_instance_of(Matrix, m.to_matrix)
		assert_instance_of(SparseMatrix, m.to_sparse_matrix)

		m = Matrix.random(4,5)
		ms = SparseMatrixCSX.from_matrix($mode, m)
		assert_equal(m, ms.to_matrix)

		m = Matrix.random(4,5)
		m[1,2] = 0
		m[3,3] = 0
		ms1 = SparseMatrix.from_matrix(m)
		ms2 = SparseMatrixCSX.from_sparse_matrix($mode, ms1)
		ms3 = ms2.to_sparse_matrix
		assert_equal(m, ms3.to_matrix)
	end

	def test_comparison
		# Should fail, not implemented
		m = SparseMatrix.new(2)
		m2 = SparseMatrix.new(2)

		assert_raise(NoMethodError) {m == m2}
		assert_raise(NoMethodError) {m =~ m2}
	end

	def test_solve
		m = Matrix[[0,1,2],[4,8,2],[2,8,6]]
		v = Vector[1,2,3]
		msc = SparseMatrixCSX.from_matrix($mode, m)

		if msc.respond_to?(:solve) && msc.type == :row

			v1 = msc.solve(v)
			v2 = (m.transpose.inverse * v).to_vector

			assert_equal(true, v1 =~ v2)
		end
	end

	def test_block_inversion
		m = Matrix.zero(9)
		m.paste!(0,0,Matrix[[10,1,2],[4,80,2],[2,8,86]])
		m.paste!(3,3,Matrix[[14,15,2],[5,7,3],[7,44,49]])
		m.paste!(6,6,Matrix[[22,4,8],[11,13,17],[8,4,62]])
		mi = m.inverse

		ms = SparseMatrixCSX.from_matrix($mode, m)
		ms.invert_diag_3x3_blocks!
		msi = ms.to_matrix

		assert_equal(true, mi =~ msi)
	end

	def test_matrix_something_operators
		# SparseMatrix * Vector
		a = Matrix[[1,2,0],[4,0,0],[0,8,9]]
		m = SparseMatrixCSX.from_matrix($mode, a)
		v = Vector[10,11,12]
		result = Vector[32.0, 40.0, 196.0]

		assert_instance_of(Vector, m.matrix_vector_multiply(v))
		assert_instance_of(Vector, m*v)

		assert_equal(result, m.matrix_vector_multiply(v))
		assert_equal(result, m*v)

		m2 = SparseMatrixCSX.from_matrix($mode, Matrix[[1,2,3],[4,5,6],[7,8,9],[10,11,12]])
		assert_equal(Vector[68.00, 167.00, 266.00, 365.00], m2*v)
		assert_raise(DimensionError) {SparseMatrixCSX.from_matrix($mode, m2.to_matrix.transpose) * v}

		m3 = SparseMatrixCSX.from_matrix($mode, Matrix[[1,2,3],[4,5,6]])
		assert_equal(Vector[68.00, 167.00], m3*v)
		assert_raise(DimensionError) {SparseMatrixCSX.from_matrix($mode, m3.to_matrix.transpose) * v}

		# SparseMatrix * Vector returning single number
		m = SparseMatrixCSX.from_matrix($mode, Matrix[[1,2,3]])
		assert_instance_of(Vector, m * m.to_matrix.transpose.to_vector)
		assert_equal(Vector[14.0], m * m.to_matrix.transpose.to_vector)
	end

	def test_in_matrix_operations
		m = Matrix[
			[0, 1, 0, 5, 1],
			[1, 2, 0, 0, 0],
			[0, 3, 0, 1, 1],
			[0, 0, 0, 5, 0],
			[0, 1, 1, 4, 0],
			[1, 2, 0, 1, 0]
		]
		ms = SparseMatrixCSX.from_matrix($mode, m)

		# column . column
		assert_equal(m.col_as_vector(1).dot(m.col_as_vector(3)), ms.col_dot_col(1,3))
		assert_equal(m.col_as_vector(0).dot(m.col_as_vector(4)), ms.col_dot_col(0,4))

		# row . row
		assert_equal(m.row_as_vector(0).dot(m.row_as_vector(4)), ms.row_dot_row(0,4))
		assert_equal(m.row_as_vector(5).dot(m.row_as_vector(1)), ms.row_dot_row(5,1))
	end

end

