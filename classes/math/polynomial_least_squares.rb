require "classes/math/polynomial.rb"
require "classes/algebra/algebra"

module PolynomialLeastSquares
	include Algebra

	def PolynomialLeastSquares.fit(x,y, order)
		k = order + 1

		# Check input data
		raise "x is not an array" unless x.kind_of?(Array)
		raise "y is not an array" unless y.kind_of?(Array)
		raise "size of x and y arrays must be the same" unless x.size == y.size

		# Matrix of the powers of x
		mx = Matrix.zero(x.size, k)

		# Matrix of y values
		my = Matrix.zero(x.size, 1)

		# Fill the matrices
		x.size.times{|i|
			k.times{|j|
				mx[i,j] = x[i]**j
			}
			my[i,0] = y[i]
		}

		# Premultiply the matrices with transpose of mx
		# to get the Vandermonde matrix
		mx_t = mx.transpose
		mx = mx_t * mx
		my = mx_t * my

		# Solve
		if order == 2
			solution = mx.inverse3 * my
		else
			begin
				solution = mx.inverse * my
			rescue
				# If the matrix elements are too small, try scaling them for more stability
				mx *= 1.0e10
				solution = mx.inverse * my
				solution *= 1.0e10
			end
		end

		return Polynomial.new(solution.column_as_array(0).reverse)
	end

end

__END__

x = [0.0, -0.007071067811865476, -0.003535533905932738, 0.003535533905932738, 0.007071067811865476]
y = [0.0, 0.06746699170734871, 0.01694649820664118, 0.015410825420940633, 0.06204432169115037]

p = PolynomialLeastSquares.fit(x,y,2)
p.print


x = [-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
y = [9.2, 4.08125, 0.85, -0.45625, 0.2, 2.85625, 7.55, 14.31875, 23.2, 34.23125, 47.45, 62.89375, 80.6]

puts PolynomialLeastSquares.fit(x,y,3).to_s
