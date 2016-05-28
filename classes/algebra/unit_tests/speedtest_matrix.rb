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


class SpeedTestMatrix < Test::Unit::TestCase

	def test_eigensystem
		# Build random, symmetric 200x200 matrix
		srand(123)
		m = Matrix.random(200)
		m = m + m.transpose
		# Diagonalize
		es = m.eigensystem
		assert_equal(3, es.size)
	end


end
