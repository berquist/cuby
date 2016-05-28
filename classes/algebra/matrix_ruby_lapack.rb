################################################################################
#
#  class Matrix
#
################################################################################

# Methods that are handled by LAPACK in the binary version

module Algebra
require @@directory+"/matrix_ruby/matrix_ruby_eigensystem"
require @@directory+"/matrix_ruby/matrix_ruby_LU"
require @@directory+"/matrix_ruby/matrix_ruby_QR"
require @@directory+"/matrix_ruby/matrix_ruby_SVD"

class Matrix
	#=======================================================================
	# determinant
	#=======================================================================
	
	def determinant
		#: Returns determinant of the matrix
		lu = LUDecomposition.new(self)
		return lu.det
	end

	#=======================================================================
	# Eigensystem
	#=======================================================================
	
	def eigensystem
		#: Solves the eigenproblem, returning matrix containing
	 	#: eigenvectors as columns (eigvs), single column matrices containing
	 	#: real (real) and imaginary (imag) parts of the eigenvalues.
	 	#| Matrix[[1,2],[3,4]].eigensystem
	 	#| => [
	 	#| | -0.82 -0.42 |
	 	#| |  0.57 -0.91 |,
	 	#| | -0.37 |
	 	#| |  5.37 |,
	 	#| |  0.00 |
	 	#| |  0.00 |]
		e = EigenvalueDecomposition.new(self)
		if e.issymmetric
			vectors = e.get_v
			real  = Matrix.zero(n,1)
			r = e.get_d
			n.times {|i| real[i,0] = r[i,i]}
			imag  = Matrix.zero(n,1)
		else
			vectors = e.get_v
			real = Matrix.columns([e.get_real_eigval])
			imag = Matrix.columns([e.get_imag_eigval])

		end

		return [vectors,real,imag]
	end

	#=======================================================================
	# LU decomposition
	#=======================================================================

	def lu
		lu = LUDecomposition.new(self)
		return [lu.get_p, lu.get_l, lu.get_u]
	end

	#=======================================================================
	# QR decomposition
	#=======================================================================
	
	def qr
		qr = QRDecomposition.new(self)
		return [qr.get_q, qr.get_r]
	end

	#=======================================================================
	# SVD
	#=======================================================================
	
	def svd
		svd = SVDecomposition.new(self)
		return[svd.get_u, svd.get_s, svd.get_v]
	end

	#=======================================================================
	# Inverse
	#=======================================================================
	
	def inverse
		#: Inverse matrix
		raise "Inverse of non-square matrices not implemented" if n != m
		### this will require to implement QR decomposition

		lu = LUDecomposition.new(self)
		return lu.solve(Matrix.identity(m))
	end

	def inverse!
		#: Inverts the matrix in place
		# Inverse
		inv = inverse
		# Copy result
		each_index{|i,j|
			self[i,j] = inv[i,j]
		}
		return nil
	end

	#=======================================================================
	# Cholesky decomposition
	#=======================================================================

	def cholesky
		#: Cholesky decomposition, returns L where M=L*L.transpose
		l = Matrix.zero(self.n)
		self.n.times{|i|
			(i+1).times{|k|
				s = 0.0
				k.times{|j| s += l[i,j] * l[k,j]; }
				if i == k
					l[i,k] = Math::sqrt(self[i,i] - s)
				else
					l[i,k] = (1.0 / l[k,k] * (self[i,k] - s))
				end
			}
		}

		return l
	end

	
end
end

