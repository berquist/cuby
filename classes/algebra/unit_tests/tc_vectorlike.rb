if ARGV[0] == 'ruby'
	$algebra_force_ruby = true
	ARGV.delete_at(0)
end

require "test/unit"
require "../algebra"
include Algebra

puts "--------------------------------------------------------------------------------"
puts "Using version: #{Algebra::version}"
puts "--------------------------------------------------------------------------------"

module Algebra
	class Veclike < Vector
		include Vector3DMethods
	end
end

class TestVector < Test::Unit::TestCase

	def test_init
		# Constructors returning Veclike
		assert_instance_of(Veclike, Veclike.from_array([1,2]))
		assert_instance_of(Veclike, Veclike.of_size(5,1.0))
		assert_instance_of(Veclike, Veclike.zero(5))

		assert_instance_of(Veclike, Veclike[1,2,3.0])
		assert_instance_of(Veclike, Veclike.from_values(1,2))
		assert_instance_of(Veclike, Veclike.random(5))

		# Zero-length should raise exception
		assert_raise(TypeError) {Veclike.from_array([])}
		assert_raise(TypeError) {Veclike.of_size(0)}
		assert_raise(TypeError) {Veclike.zero(0)}
		assert_raise(TypeError) {Veclike[]}
		assert_raise(TypeError) {Veclike.from_values()}
		assert_raise(TypeError) {Veclike.random(0)}

		# Negative length should raise exception
		assert_raise(TypeError) {Veclike.of_size(-1)}
		assert_raise(TypeError) {Veclike.zero(-1)}
		assert_raise(TypeError) {Veclike.random(-1)}

		# Size and value check
		assert_equal(Veclike.of_size(5,1.0).size, 5)
		assert_equal(Veclike.of_size(5,1.0)[3], 1.0)
	end

	def test_basic
		# size
		assert_equal(Veclike[1,2,3.0].size,3)

		# copy_from!
		v1 = Veclike[1,2,3]
		v2 = Veclike[4,5,6]
		v3 = Veclike[1,5]
		v1_id = v1.object_id
		v1.copy_from!(v2)
		assert_equal(v1,v2)
		assert_equal(v1_id,v1.object_id)
		v1.copy_from!([4,5,6])
		assert_equal(v1,v2)
		assert_raise(DimensionError) {v1.copy_from!(v3)}

		# clone / deep_copy
		v1 = Veclike[1,2,3]
		v2 = v1.clone
		assert_equal(v1,v2)
		assert_not_equal(v1.object_id,v2.object_id)
		v2[0] = 5
		assert_not_equal(v1,v2)
		v2 = v1.deep_copy
		assert_equal(v1,v2)
	end

	def test_element_access
		# Read
		v1 = Veclike[1,2,3]
		assert_equal(v1[1],2.0)
		assert_raise(IndexError) {v1[3]}
		assert_raise(IndexError) {v1[-1]}

		# Write
		v1[2] = 5.0
		assert_equal(v1[2],5.0)
		assert_raise(IndexError) {v1[3] = 0.0}
		assert_raise(IndexError) {v1[-1] = 0.0}
		assert_raise(TypeError) {v1[2] = "foo"}

		# Subvector
		v1 = Veclike[1,2,3,4,5,6,7,8,9]
		assert_equal(v1.subvector(4,3), Veclike[5,6,7])
		assert_raise(IndexError) {v1.subvector(7,3)}
		assert_raise(IndexError) {v1.subvector(-1,3)}

		# first, last
		assert_equal(1.0, v1.first)
		assert_equal(9.0, v1.last)
	end

	def test_type_conversion
		# to_a
		v1 = Veclike[1,2,3]
		assert_instance_of(Array, v1.to_a)
		# to_matrix
		assert_instance_of(Matrix, v1.to_matrix)
		assert_equal(v1.to_matrix.m, 3)
		assert_equal(v1.to_matrix.n, 1)
		# to_s, inspect
		assert_instance_of(String, v1.to_s)
		assert_instance_of(String, v1.inspect)
		# to_vector
		assert_instance_of(Veclike, v1.to_vector)
		assert_not_equal(v1.object_id,v1.to_vector.object_id)
	end

	def test_iterators
		v1 = Veclike[1,2,3.33]
		# each_index
		count = 0
		last = nil
		v1.each_index { |i| count += 1; last = i}
		assert_equal(count, 3)
		assert_equal(last, 2)
		
		#each
		count = 0
		last = nil
		v1.each { |i| count += 1; last = i}
		assert_equal(count, 3)
		assert_equal(last, 3.33)

		# each_with_index
		vals = []
		indices = []
		v1.each_with_index { |v,i|
			vals << v
			indices << i
		}
		assert_equal(vals, [1,2,3.33])
		assert_equal(indices, [0,1,2])

		# collect
		result = v1.collect{|i| i.to_s}
		assert_equal(result,["1.0","2.0","3.33"])
	end

	def test_unary_operators
		v1 = Veclike[1,2,3]
		assert_equal(v1,+v1)
		assert_not_equal(v1,-v1)
		assert_equal(v1,-v1 * -1.0)
		assert_equal(Veclike[-1,-2,-3],-v1)

		assert_instance_of(Veclike, +v1)
		assert_instance_of(Veclike, -v1)
	end

	def test_vector_vector_operators
		# + and -
		assert_equal(Veclike[1,2] + Veclike[3,4], Veclike[4,6])
		assert_instance_of(Veclike, Veclike[1,2] + Veclike[3,4])
		assert_raise(DimensionError) {Veclike[1,2] + Veclike[3,4,5]}

		assert_equal(Veclike[1,2] - Veclike[0,1], Veclike[1,1])
		assert_instance_of(Veclike, Veclike[1,2] - Veclike[3,4])
		assert_raise(DimensionError) {Veclike[1,2] - Veclike[3,4,5]}

		# plus!	
		v1 = Veclike[1,2]
		old_id = v1.object_id
		v1.plus!(Veclike[3,4])	
		assert_equal(v1, Veclike[4,6])
		assert_instance_of(Veclike, v1)
		assert_equal(old_id, v1.object_id)

		# minus!
		v1 = Veclike[1,2]
		old_id = v1.object_id
		v1.minus!(Veclike[0,1])	
		assert_equal(v1, Veclike[1,1])
		assert_instance_of(Veclike, v1)
		assert_equal(old_id, v1.object_id)

		# dot_product (synonym dot)
		assert_equal(Veclike[1,0].dot_product(Veclike[0,1]), 0.0)
		assert_equal(Veclike[1,0].dot(Veclike[0,1]), 0.0)
		assert_equal(Veclike[1,1].dot(Veclike[2,1]), 3.0)
		assert_instance_of(Float,Veclike[1,0].dot(Veclike[0,1]))

		# elementwise multiply and divide
		v = Veclike[1,2,3]
		assert_equal(Veclike[2,3,9], v.elementwise_multiply(Veclike[2,1.5,3]))
		assert_equal(Veclike[0.5,4,1], v.elementwise_divide(Veclike[2,0.5,3]))
	end

	def test_vector_something_operators
		# multiplication
		assert_instance_of(Float, Veclike[1,0] * Veclike[0,1])
		assert_instance_of(Veclike, Veclike[1,0] * 2)
		assert_instance_of(Veclike, Veclike[1,0] * 2.0)
		assert_raise(TypeError) {Veclike[1,0] * "b"}
		assert_equal( Veclike[1,0] * Veclike[0,1], 0.0 )
		assert_equal(Veclike[1,0.5] * 2, Veclike[2,1])
		assert_equal(Veclike[1,0.5] * 2.0, Veclike[2,1])

		#division
		assert_instance_of(Veclike, Veclike[1,0] / 2)
		assert_instance_of(Veclike, Veclike[1,0] / 2.0)
		assert_raise(TypeError) {Veclike[1,0] / "b"}
		assert_raise(TypeError) {Veclike[1,0] / Veclike[2,1]}
		assert_equal(Veclike[4,2] / 2, Veclike[2,1])
		assert_equal(Veclike[4,2] / 2.0, Veclike[2,1])
	end

	def test_comparison
		assert_equal(false, Veclike[1,0] == "a")
		assert_raise(DimensionError) {Veclike[1,0] == Veclike[1,0,2]}
		assert_equal(true, Veclike[1,0] == Veclike[1,0])
		assert_equal(false, Veclike[1,0] == Veclike[1,1e-30])
		assert_equal(false, Veclike[1,0] == Veclike[1,5])

		Veclike.epsilon = 1.0e-8
		assert_equal(Veclike[0,1] =~ Veclike[0.0,1.0+1.0e-9], true)
		assert_equal(Veclike[0,1] =~ Veclike[0.0,1.0+1.0e-7], false)

	end

	def test_misc
		# absolute, abs, r, rms
		assert_equal(5.0, Veclike[3,4].absolute)
		assert_equal(5.0, Veclike[3,4].abs)
		assert_equal(5.0, Veclike[3,4].r)
		assert_equal(5.0/2**0.5, Veclike[3,4].rms)


		# normalize
		assert_equal(Veclike[0,4].normalize, Veclike[0, 1])
		assert_in_delta(2.0**0.5/2, Veclike[1,1,0,0,0].normalize[0], 1e-10)
		assert_instance_of(Veclike, Veclike[2,2].normalize)

		# normalize!
		v = Veclike[0,2]
		old_id = v.object_id
		v.normalize!
		assert_equal(v, Veclike[0,1])
		assert_equal(v.object_id, old_id)

		# distance, dist
		assert_in_delta(2.0**0.5, Veclike[1,0].distance(Veclike[2,1]), 1e-20)
		assert_in_delta(2.0**0.5, Veclike[1,0].dist(Veclike[2,1]), 1e-20)

		# sum
		v = Veclike[0,1,3,-6,-6.1,5]
		assert_in_delta(-3.1, v.sum, 1.0e-10)
		# max, min
		assert_equal(-6.1, v.min)
		assert_equal(5.0, v.max)
		# max_abs, min_abs
		assert_equal(0.0, v.min_abs)
		assert_equal(6.1, v.max_abs)
	end

	def test_yaml
		# Psych version
		s1 = "--- !cuby.molecular.cz,2009/Algebra::Veclike {elements: [1.0, 2.0, 3.0, 4.0, 5.0]}\n"
		# Syck version
		s2 = "--- !cuby.molecular.cz,2009/Algebra::Veclike \nelements: [1.0, 2.0, 3.0, 4.0, 5.0]\n"
		v = Veclike[1,2,3,4,5]

		assert_equal(v, YAML.load(s1))
		assert_equal(v, YAML.load(s2))

		assert_equal(v, YAML.load(v.to_yaml))
	end

	def test_3D_only
		# cross product
		assert_equal(Veclike[1,0,0].cross_product(Veclike[0,1,0]), Veclike[0.0, 0.0, 1.0])
		assert_instance_of(Veclike, Veclike[1,0,0].cross_product(Veclike[0,1,0]))
		assert_raise(Not3DError) {Veclike[1,0,0].cross_product(Veclike[0,1])}
		assert_raise(Not3DError) {Veclike[1,0].cross_product(Veclike[0,1,0])}

		# angle
		assert_in_delta(Math::PI/2, Veclike[1,0,0].angle(Veclike[0,1,0]), 1e-10)
		assert_equal(Veclike[1,0,0].angle(Veclike[1,0,0]), 0.0)

		# angle_deg
		assert_in_delta(90.0, Veclike[1,0,0].angle_deg(Veclike[0,1,0]), 1e-10)

		# rotate_angle_axis
		# (tests also the internally used construction of rotator matrix and rotation using it)
		v = Veclike[1,2,3]
		v2 = v.rotate_angle_axis(Math::PI / 180 * 30, Veclike[0.5,-1,0.4])
		assert_in_delta(v.abs, v2.abs, 1e-10)
		assert_in_delta(-0.748315862860624, v2.x, 1e-10)
		assert_in_delta(1.29737248421965, v2.y, 1e-10)
		assert_in_delta(3.42882603912489, v2.z, 1e-10)
		assert_instance_of(Veclike, v2)
		assert_raise(Not3DError) {v.rotate_angle_axis(Math::PI / 180 * 30, Veclike[0.5,-1])}
		assert_raise(Not3DError) {Veclike[1,2].rotate_angle_axis(Math::PI / 180 * 30, Veclike[0.5,-1,0.4])}

		# rotate_angle_axis!
		# (tests also the internally used construction of rotator matrix and the in-place rotation using it)
		v = Veclike[1,2,3]
		v.rotate_angle_axis!(Math::PI / 180 * 30, Veclike[0.5,-1,0.4])
		assert_in_delta(-0.748315862860624, v.x, 1e-10)
		assert_in_delta(1.29737248421965, v.y, 1e-10)
		assert_in_delta(3.42882603912489, v.z, 1e-10)
		assert_instance_of(Veclike, v)

		# direct access to coordinates
		assert_raise(Not3DError) {Veclike[1,2,3,4].x}
		assert_raise(Not3DError) {Veclike[1,2,3,4].x=5}
		v = Veclike[1,2,3]
		v.x = 4
		v.y = 5
		v.z = 6
		assert_equal(v.x, 4)
		assert_equal(v.y, 5)
		assert_equal(v.z, 6)
	end

	def test_float_extension
		assert_instance_of(Veclike, 2.0 * Veclike[1,0])
		assert_equal(2.0 * Veclike[1,0.5], Veclike[2,1])

	end

	def test_fixnum_extension
		assert_instance_of(Veclike, 2 * Veclike[1,0])
		assert_equal(2 * Veclike[1,0.5], Veclike[2,1])
	end
end
