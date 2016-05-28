################################################################################
#
#  class Matrix : Fixes to C version
#
################################################################################

module Algebra
class Matrix
	def qr
		#: QR decomposition, returns [Q,R]
		x,r = self.qr2
		q = Matrix.diagonal(x.m)
		i = Matrix.diagonal(x.m)

		[m,n].min.times{|j|
			v = x.submatrix(0,j,x.m,1)
			q = q * (i - v * v.transpose)
		}
		return [q,r]
	end
end
end
