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

module Algebra
	class Matlike < Matrix
	end
end

class TestMatlike < Test::Unit::TestCase

	def test_init
		# Constructors, should return instance of Matlike
		assert_instance_of(Matlike, Matlike[[1,2],[3,4]])
		assert_instance_of(Matlike, Matlike.from_array([[1,2],[3,4]]))
		assert_instance_of(Matlike, Matlike.rows([[1,2],[3,4]]))
		assert_instance_of(Matlike, Matlike.columns([[1,3],[2,4]]))
		assert_instance_of(Matlike, Matlike.zero(5))
		assert_instance_of(Matlike, Matlike.zero(5,10))
		assert_instance_of(Matlike, Matlike.diagonal(5))
		assert_instance_of(Matlike, Matlike.diagonal(5,5))
		assert_instance_of(Matlike, Matlike.diagonal(5,6))
		assert_instance_of(Matlike, Matlike.diagonal([1,2,3]))
		assert_instance_of(Matlike, Matlike.diagonal(Vector[1,2,3]))
		assert_instance_of(Matlike, Matlike.identity(5))
		assert_instance_of(Matlike, Matlike.filled(4,5,1.00))
		assert_instance_of(Matlike, Matlike.random(2,3))

		# Check size of the matrix created
		m = Matlike[[1,2],[3,4],[5,6]]
		assert_equal(3,m.m)
		assert_equal(2,m.n)
		m = Matlike.identity(5)
		assert_equal(5,m.m)
		assert_equal(5,m.n)

		# Check contents of some matrices
		assert_equal(Matlike.identity(5), Matlike.diagonal(5,1))
		assert_equal(Matlike.rows([[1,2],[3,4]]), Matlike.columns([[1,3],[2,4]]))

		# Check error handling in Matlike.from_array
		assert_raise(TypeError) { Matlike.from_array([])}
		assert_raise(TypeError) { Matlike.from_array([[],[]])}
		assert_raise(TypeError) { Matlike.from_array("Foo")}
		assert_raise(TypeError) { Matlike.from_array([[1,2],"Foo"])}
		assert_raise(TypeError) { Matlike.from_array([[1,2],["Foo",4]])}
		assert_raise(TypeError) { Matlike.from_array([[1,2],[3,4,5]])}

		# Zero-sized matrices
		assert_raise(TypeError) { Matlike.zero(0,0) }
		assert_raise(TypeError) { Matlike.zero(1,0) }
		assert_raise(TypeError) { Matlike.zero(0,1) }
		assert_raise(TypeError) { Matlike.diagonal(0,1) }
		assert_raise(TypeError) { Matlike.diagonal([]) }
		assert_raise(TypeError) { Matlike.identity(0) }
		assert_raise(TypeError) { Matlike.filled(0,0,1) }
		assert_raise(TypeError) { Matlike.filled(1,0,1) }
		assert_raise(TypeError) { Matlike.filled(0,1,1) }
		assert_raise(TypeError) { Matlike.random(0,0) }
		assert_raise(TypeError) { Matlike.random(1,0) }
		assert_raise(TypeError) { Matlike.random(0,1) }

		# Negative-sized matrices
		assert_raise(TypeError) { Matlike.zero(-1,0) }
		assert_raise(TypeError) { Matlike.zero(0,-1) }
		assert_raise(TypeError) { Matlike.diagonal(-1,1) }
		assert_raise(TypeError) { Matlike.identity(-1) }
		assert_raise(TypeError) { Matlike.filled(-1,0,1) }
		assert_raise(TypeError) { Matlike.filled(0,-1,1) }
		assert_raise(TypeError) { Matlike.random(-1,0) }
		assert_raise(TypeError) { Matlike.random(0,-1) }
	end

	def test_basic
		# Dimensions
		m = Matlike[[1,2],[3,4],[5,6]]
		assert_equal(3,m.m)
		assert_equal(2,m.n)
		assert_equal([3,2],m.dimensions)
		assert_equal([3,2],m.size)
		assert_equal("3 x 2",m.size_str)

		# Clone / deep_copy
		m = Matlike[[1,2],[3,4]]
		m2 = m.clone
		assert_equal(m,m2)
		assert_not_equal(m.object_id,m2.object_id)
		m2[1,1] = 0
		assert_not_equal(m,m2)
	end

	def test_element_access
		# Read/write element
		m = Matlike[[1,2],[3,4]]
		assert_equal(1,m[0,0])
		assert_equal(2,m[0,1])
		assert_equal(3,m[1,0])
		assert_equal(4,m[1,1])
		m[1,0] = 15
		assert_equal(15,m[1,0])

		# Exceptions on read/write
		assert_raise(IndexError) {m[0,3]}
		assert_raise(IndexError) {m[3,0]}
		assert_raise(IndexError) {m[-0,3]}
		assert_raise(IndexError) {m[-3,0]}
		assert_raise(IndexError) {m[0,3] = 15}
		assert_raise(IndexError) {m[3,0] = 12}
		assert_raise(IndexError) {m[0,-3] = 15}
		assert_raise(IndexError) {m[-3,0] = 12}
		assert_raise(TypeError)  {m[0,0] = "foo"}
	end

	def test_row_col_submatrix_access
		# Row/col access
		m = Matlike[[2,4,6],[3,5,7]]
		assert_instance_of(Array, m.row_as_array(1))
		assert_equal([3.0,5.0,7.0], m.row_as_array(1))

		assert_instance_of(Array, m.col_as_array(1))
		assert_equal([4.0,5.0], m.col_as_array(1))
		assert_equal([4.0,5.0], m.column_as_array(1))

		assert_instance_of(Vector, m.row_as_vector(1))
		assert_equal(Vector[3.0,5.0,7.0], m.row_as_vector(1))

		assert_instance_of(Vector, m.col_as_vector(1))
		assert_equal(Vector[4.0,5.0], m.col_as_vector(1))
		assert_equal(Vector[4.0,5.0], m.column_as_vector(1))

		# Submatrix
		m = Matlike[[1.0, 2.0, 3.0, 4.0, 5.0],
		           [6.0, 7.0, 8.0, 9.0, 8.0],
			   [7.0, 6.0, 5.0, 4.0, 3.0],
			   [2.0, 1.0, 0.0, 1.0, 2.0]]

		sm = m.submatrix(1,2,2,3)
		assert_instance_of(Matlike, sm)
		assert_equal(2, sm.m)
		assert_equal(3, sm.n)
		assert_equal(Matlike[[8.0, 9.0, 8.0],[5.0, 4.0, 3.0]], sm)
		assert_raise(IndexError) {m.submatrix(1,2,2,9)}
		assert_raise(IndexError) {m.submatrix(1,2,9,2)}
		assert_raise(IndexError) {m.submatrix(9,2,1,2)}
		assert_raise(IndexError) {m.submatrix(1,9,1,2)}

		# Paste!
		mp = Matlike[[0.1, 0.2],[0.3, 0.4]]
		m.paste!(1,3,mp)
		assert_equal(0.1, m[1,3])
		assert_equal(0.4, m[2,4])
		assert_raise(IndexError) {m.paste!(3,1,mp)}
		assert_raise(IndexError) {m.paste!(0,4,mp)}

		# Swap_columns!
		m = Matlike[[1,2,3,4,5,6],[1.1,2.2,3.3,4.4,5.5,6.6]]
		msw = Matlike[[6,2,3,4,5,1],[6.6,2.2,3.3,4.4,5.5,1.1]]

		msw.swap_columns!(0,5)
		assert_equal(m, msw)

		# Reorder_columns!
		mr = Matlike[[1,2,6,4,5,3]]
		mr.reorder_columns!([0,1,5,3,4,2])
		assert_equal([1.0,2.0,3.0,4.0,5.0,6.0], mr.row_as_array(0))

	end

	def test_type_conversion
		m = Matlike[[1,2],[3,4]]
		assert_instance_of(String, m.to_s)
		assert_instance_of(String, m.inspect)
		assert_instance_of(String, m.to_s_rb)

		assert_instance_of(Array, m.to_a)
		assert_equal([[1.0,2.0],[3.0,4.0]], m.to_a)
		assert_instance_of(Array, m.to_a[0])

		assert_raise(DimensionError) {m.to_vector}
		m2 = Matlike[[1],[2],[3]]
		assert_instance_of(Vector, m2.to_vector)
		assert_equal(Vector[1,2,3], m2.to_vector)
	end

	def test_iterators
		m = Matlike[[1,2],[3,4]]
		# each
		s = ''
		m.each {|e| s += e.to_s + " "}
		assert_equal("1.0 2.0 3.0 4.0 ", s)
		# each_index
		s = ''
		m.each_index {|i,j| s += "#{i},#{j} "}
		assert_equal("0,0 1,0 0,1 1,1 ", s)
		#each_with_index
		s = ''
		m.each_with_index {|e,i,j| s += "#{i},#{j}=#{e} "}
		assert_equal("0,0=1.0 1,0=3.0 0,1=2.0 1,1=4.0 ", s)
		#row_each_with_index
		av = []
		aj = []
		m.row_each_with_index(0){|e,j|
			av << e
			aj << j
		}
		assert_equal([1,2], av)
		assert_equal([0,1], aj)
		# each_diagonal
		a = []
		m.each_diagonal{|x| a << x}
		assert_equal([1,4],a)
	end

	def test_min_max_norms
		m = Matlike[[1,-2],[3,4],[0.2,1.5]]
		assert_equal(-2.0, m.min)
		assert_equal(4.0, m.max)
		assert_equal([-2.0,0,1], m.min_with_index)
		assert_equal([4.0,1,1], m.max_with_index)
	end

	def test_unary_operators
		m = Matlike[[1,2],[3,4]]
		assert_equal(m,+m)
		assert_not_equal(m,-m)
		assert_equal(m,-m * -1.0)

		assert_instance_of(m.class, +m)
		assert_instance_of(m.class, -m)
	end

	def test_matrix_matrix_operators
		m1 = Matlike[[1,2],[3,4]]
		m2 = Matlike[[1,0],[0,2]]
		# +
		assert_instance_of(Matlike, m1 + m2)
		assert_equal(Matlike[[2,2],[3,6]], m1 + m2)
		assert_raise(DimensionError) {m1 + Matlike[[1,0]]}
		assert_raise(TypeError) {m1 + 1}
		# -
		assert_instance_of(Matlike, m1 + m2)
		assert_equal(Matlike[[0,2],[3,2]], m1 - m2)
		assert_raise(DimensionError) {m1 - Matlike[[1,0]]}
		# plus!
		m1 = Matlike[[1,2],[3,4]]
		m1_id = m1.object_id
		m1.plus! Matlike[[1,0],[0,2]]
		assert_equal(Matlike[[2,2],[3,6]], m1)
		assert_equal(m1_id, m1.object_id)
		assert_raise(DimensionError) {m1.plus! Matlike[[1,0]]}
		# minus!
		m1 = Matlike[[1,2],[3,4]]
		m1_id = m1.object_id
		m1.minus! Matlike[[1,0],[0,2]]
		assert_equal(Matlike[[0,2],[3,2]], m1)
		assert_equal(m1_id, m1.object_id)
		assert_raise(DimensionError) {m1.minus! Matlike[[1,0]]}
		# mul!
		m1 = Matlike[[1,2],[3,4]]
		m1_id = m1.object_id
		m1.mul!(2)
		assert_equal(Matlike[[2,4],[6,8]], m1)
		assert_equal(m1_id, m1.object_id)
		assert_raise(TypeError) {m1.mul!(Matlike[[1,2],[3,4]])}
		# div!
		m1 = Matlike[[1,2],[3,4]]
		m1_id = m1.object_id
		m1.div!(2)
		assert_equal(Matlike[[0.5,1],[1.5,2]], m1)
		assert_equal(m1_id, m1.object_id)
		# m * m (all size combinations)
		assert_equal(Matlike[[1,2],[3,4],[5,6],[7,8]] * Matlike[[1,2,3],[4,5,6]],
			   Matlike[[  9.00, 12.00, 15.00 ],
				  [ 19.00, 26.00, 33.00 ],
				  [ 29.00, 40.00, 51.00 ],
				  [ 39.00, 54.00, 69.00 ]])
		assert_equal(Matlike[[1,2],[3,4],[5,6]] * Matlike[[1,2,3],[4,5,6]],
			     Matlike[[ 9.00, 12.00, 15.00],
			            [19.00, 26.00, 33.00],
			            [29.00, 40.00, 51.00]])
		assert_equal(Matlike[[1,2],[3,4]] * Matlike[[1,2,3],[4,5,6]],
			     Matlike[[ 9.00, 12.00, 15.00],
			            [19.00, 26.00, 33.00]])
		assert_equal(Matlike[[1,2],[3,4]] * Matlike[[1,2],[7.1,8.5]],
			     Matlike[[15.20, 19.00],
			            [31.40, 40.00]])
		assert_raise(DimensionError) {
			Matlike[[1,2],[3,4],[5,6]] * Matlike[[1,2],[3,4],[4,5]]
		}

		# m * m returning single number
		m = Matlike[[1,2,3]]
		assert_instance_of(Matlike, m * m.transpose)
		assert_equal(Matlike[[14.0]], m * m.transpose)
	end

	def test_matrix_something_operators
		# Matlike * Vector
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		v = Vector[10,11,12]
		assert_instance_of(Vector, m.matrix_vector_multiply(v))
		assert_instance_of(Vector, m*v)
		assert_equal(Vector[68.00, 167.00, 266.00], m.matrix_vector_multiply(v))
		assert_equal(Vector[68.00, 167.00, 266.00], m*v)

		m2 = Matlike[[1,2,3],[4,5,6],[7,8,9],[10,11,12]]
		assert_equal(Vector[68.00, 167.00, 266.00, 365.00], m2*v)
		assert_raise(DimensionError) {m2.transpose * v}

		m3 = Matlike[[1,2,3],[4,5,6]]
		assert_equal(Vector[68.00, 167.00], m3*v)
		assert_raise(DimensionError) {m3.transpose * v}

		# Matlike * Vector returning single number
		m = Matlike[[1,2,3]]
		assert_instance_of(Vector, m * m.transpose.to_vector)
		assert_equal(Vector[14.0], m * m.transpose.to_vector)

		# Matlike * Numeric
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		m2 = Matlike[[ 2.00,  4.00,  6.00], [ 8.00, 10.00, 12.00], [14.00, 16.00, 18.00]]

		assert_equal(m2, m*2.0)
		assert_equal(m2, m*2)
		assert_raise(TypeError) {m * "x"}

		# Matlike / Numeric
		assert_equal(m, m2/2)
		assert_equal(m, m2/2.0)
		assert_raise(TypeError) {m / m2}
		assert_raise(TypeError) {m / v}

		# vec.T * mat * vec
		# returns number
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		v1 = Vector[2,5,8]
		v2 = Vector[9,6,3]
		assert_equal(m.vt_self_v(v1,v2), 1584.0)
	end

	def test_float_extension
		assert_instance_of(Matlike, 2.0 * Matlike[[1,0]])
		assert_equal(Matlike[[2,1]], 2.0 * Matlike[[1,0.5]])

	end

	def test_fixnum_extension
		assert_instance_of(Matlike, 2 * Matlike[[1,0]])
		assert_equal(2 * Matlike[[1,0.5]], Matlike[[2,1]])
	end

	def test_comparison
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		m2 = m.deep_copy
		m3 = Matlike[[1,2,3],[4,5,6]]
		mx = Matlike[[1,2,3],[4,5,6],[8,8,9]]

		assert_equal(true, m == m2)
		assert_equal(false, m == mx)
		assert_raise(DimensionError) {m == m3}

		# Loose comparison
		m2[0,0] += Matlike::epsilon / 2.0
		assert_equal(true, m =~ m2)
		assert_equal(false, m == m2)
		Matlike::epsilon /= 3
		assert_equal(false, m =~ m2)
		Matlike::epsilon = Matlike::DEFAULT_EPSILON # restore
	end

	def test_matrix_algebra
		# Transpose
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		mt = Matlike[[1,4,7],[2,5,8],[3,6,9]]
		assert_equal(m, mt.transpose)
		mt.transpose!
		assert_equal(m, mt)

		# Determinant
		assert_equal(-2.0, Matlike[[1,2],[3,4]].det)
		assert_equal(-2.0, Matlike[[1,2],[3,4]].determinant)
		assert_equal(300.0, Matlike[[-1,2,3],[3,-5,6],[7,8,9]].det)
		assert_raise(DimensionError) {Matlike[[-1,2,3],[3,-5,6]].det}

		# Eigenvalues
		m = Matlike[[-1, 3, -1],[-3,5,-1],[-3,3,1]]
		vec, real, imag = m.eigensystem
		assert_equal([3,3], vec.size)
		assert_equal([3,1], real.size)
		assert_equal([3,1], imag.size)
		assert_equal(0.0, imag.max)
		diagev = Matlike.diagonal(real.to_vector)
		assert_equal(true, (vec.inverse * m * vec) =~ diagev)

		# Eigenvalue sorter
		m = Matlike[[ 0.41,  0.16,  0.63,  0.91,  0.01],
			   [ 0.00,  0.32,  0.96,  0.52,  0.82],
			   [ 0.53,  0.07,  0.01,  0.56,  0.11],
			   [ 0.53,  0.04,  0.47,  0.28,  0.96],
			   [ 0.69,  0.61,  0.90,  0.67,  0.40]]
		m = m + m.transpose
		[:ascending_real, :descending_real, :ascending_real_abs, :descending_real_abs].each{|order|
			vec, real, imag = m.eigensystem
			Matlike.sort_eigensystem!(vec, real, imag, order)
			diagev = Matlike.diagonal(real.to_vector)
		        assert_equal(true, (vec.inverse * m * vec) =~ diagev)
		}

		# LU
		# square
		m = Matlike[[1,2,3],[0,1,4],[5,6,0]]
		p,l,u = m.lu
		assert_equal(true, m =~ p*l*u)
		# 4 * 3
		m = Matlike[[ 0.62,  0.50,  0.68], [ 0.53,  0.06,  0.51], [ 0.83,  0.78,  0.28], [ 0.06,  0.08,  0.25]]
		p,l,u = m.lu
		assert_equal(true, m =~ p*l*u)
		# 3 * 4
		m = m.transpose
		p,l,u = m.lu
		assert_equal(true, m =~ p*l*u)
		
		# QR
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		q,r = m.qr
		assert_equal(true, m =~ q*r)

		# SVD
		u,s,v = m.svd
		assert_equal(true, u * s * v =~ m)

		# Inverse
		m2 = Matlike[[2,3],[2,2]]
		assert_equal(true, m2 * m2.inverse =~ Matlike.unit(2))
		assert_equal(true, m2 * m2.inverse2 =~ Matlike.unit(2))

		m3 = Matlike[[1,2,3],[0,1,4],[5,6,0]]
		assert_equal(true, m3 * m3.inv =~ Matlike.unit(3))
		assert_equal(true, m3 * m3.inverse3 =~ Matlike.unit(3))

		m4 = m3.deep_copy
		m4.inv!
		assert_equal(true, m3 * m4 =~ Matlike.unit(3))

		# Inverse: numerical stability test for nearly-singular matrices
		assert_raise(RuntimeError) {Matlike[[1,2,3],[4,5,6],[7,8,9.0000001]].inverse}
		assert_nothing_raised(RuntimeError) {Matlike[[1,2,3],[4,5,6],[7,8,9.000001]].inverse}

		# Moore-Penrose pseudoinverse
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		assert_equal(true, m * m.pseudoinverse * m =~ m)
		m = Matlike[[1,2,3],[4,5,6],[7,8,9],[1,2,3]]
		assert_equal(true, m * m.pseudoinverse * m =~ m)
		m = Matlike[[1,2,3],[4,5,6]]
		assert_equal(true, m * m.pseudoinverse * m =~ m)

		# Trace
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		assert_equal(15, m.trace)

	end

	def test_yaml
		# Read own output
		m = Matlike[[1,2,3],[4,5,6],[7,8,9]]
		s = m.to_yaml
		m2 = YAML.load(s)
		assert_equal(true, m =~ m2)

		# Read external output
		# from ruby1.9
		m2 = YAML.load("--- !cuby.molecular.cz,2009/Algebra::Matlike\nelements:\n- - 1.0\n  - 2.0\n  - 3.0\n- - 4.0\n  - 5.0\n  - 6.0\n- - 7.0\n  - 8.0\n  - 9.0\n")
		assert_equal(true, m =~ m2)
		# from ruby1.8
		m2 = YAML.load("--- !cuby.molecular.cz,2009/Algebra::Matlike \nelements: \n- [1.0, 2.0, 3.0]\n- [4.0, 5.0, 6.0]\n- [7.0, 8.0, 9.0]\n")
		assert_equal(true, m =~ m2)
	end

	def test_input_output
		# File read/write
		m = Matlike.random(3)
		m.to_file("test_mat.txt")
		m2 = Matlike.from_file("test_mat.txt")
		FileUtils.rm("test_mat.txt")
		assert_equal(true, m =~ m2)
	end

end
