module Algebra
class LUDecomposition 

	def initialize(a) # a : Matrix
		@lu = a.to_a # array of arrays
		@m = a.m
		@n = a.n
		@piv = []; @m.times {|i| @piv << i}
		@pivsign = 1.0
		lurowi = []
		lucolj = []

		# outer loop
		@n.times{|j|
			# copy j-th column
			@m.times{|i| lucolj[i] = @lu[i][j]}

			# apply previous transformation
			@m.times{|i|
				lurowi = @lu[i]
				kmax = [i,j].min
				s = 0.0
				kmax.times{|k|
					s += lurowi[k] * lucolj[k]
				}
				lurowi[j] = lucolj[i] -= s
			}

         		# Find pivot and exchange if necessary
			p = j
			(@m - j - 1).times{|ii|	# for (int i = j+1; i < m; i++)
				i = ii + j + 1
				if lucolj[i].abs > (lucolj[p])
					p = i
				end
			}
			if p != j
				@n.times {|k|
					t = @lu[p][k]; @lu[p][k] = @lu[j][k]; @lu[j][k] = t
				}
				k = @piv[p]; @piv[p] = @piv[j]; @piv[j] = k
				@pivsign = -@pivsign
			end

			# Compute multipliers
			if j < @m && @lu[j][j] != 0.0
				(@m - j - 1).times{|ii|
					i = ii + j + 1
					@lu[i][j] /= @lu[j][j]
				}
			end
		}
	end

	def is_nonsingular? # true if U, and hence A, is nonsingular
		@n.times{|j|
			return false if @lu[j][j] == 0.0
		}
		return true
	end

	def raise_if_nearly_singular
		@n.times{|j|
			raise "Matrix is nearly singular, solution numerically unstable" if @lu[j][j].abs < 1.0e-7
		}
		return nil
	end

	def get_l # lower triangular factor
		minmn = [@m,@n].min
		l = Matrix.zero(@m,minmn)
		@m.times{|i|
			minmn.times{|j|
				if i > j
					l[i,j] = @lu[i][j]
				elsif i == j
					 l[i,j] = 1.0
				end
			}
		}
		return l
	end

	def get_u
		minmn = [@m,@n].min
		u = Matrix.zero(minmn,@n)
		minmn.times{|i|
			@n.times{|j|
				if i <= j
					u[i,j] = @lu[i][j]
				end
			}
		}
		return u
	end

	def get_p
		p = Matrix.zero(@m,@m)
		@m.times {|i|
			p[@piv[i],i] = 1.0
		}
		return p
	end

	def det
		raise(DimensionError,"Determinant can be calculated only for square matrix") if @m != @n
		d = @pivsign
		@n.times{|j|
			d *= @lu[j][j]
		}
		return d
	end

	def solve(b) # => Matrix
		raise(DimensionError,"Number of rows must be identical")if b.m != @m
		raise "Matrix is singular" unless is_nonsingular?
		raise_if_nearly_singular

		# Copy right hand side with pivoting
		nx = b.m # B.getColumnDimension();
		# matrix with swapped rows made from b
		#Matrix Xmat = B.getMatrix(piv,0,nx-1);
		#double[][] X = Xmat.getArray();
		xmat = Matrix.zero(nx, @piv.size)
		row = 0
		@piv.each{|p|
			p = p.to_i
			nx.times {|col|
				xmat[row,col] = b[p,col]
			}
			row += 1
		}

		# Solve L*Y = B(piv,:)
		@n.times{|k| # for (int k = 0; k < n; k++)
			(@n-k-1).times {|ii| # for (int i = k+1; i < n; i++)
				i = ii + k + 1
				nx.times{|j| # for (int j = 0; j < nx; j++)
					# X[i][j] -= X[k][j]*LU[i][k];
					xmat[i,j] -= xmat[k,j] * @lu[i][k]
				}
			}
		}
		# Solve U*X = Y;
		@n.times {|kk| # for (int k = n-1; k >= 0; k--)
			k = @n-1-kk
			nx.times{|j| # for (int j = 0; j < nx; j++)
				# X[k][j] /= LU[k][k];
				xmat[k,j] /= @lu[k][k]
			}
			k.times {|i| #for (int i = 0; i < k; i++)
				nx.times{|j| # for (int j = 0; j < nx; j++)
					#X[i][j] -= X[k][j]*LU[i][k];
					xmat[i,j] -= xmat[k,j] * @lu[i][k]
				}
			}
		}
		return xmat
	end

end
end


#--------------------------------------------------
# 
#    /** Determinant
#    @return     det(A)
#    @exception  IllegalArgumentException  Matrix must be square
#    */
# 
#    public double det () {
#       if (m != n) {
#          throw new IllegalArgumentException("Matrix must be square.");
#       }
#       double d = (double) pivsign;
#       for (int j = 0; j < n; j++) {
#          d *= LU[j][j];
#       }
#       return d;
#    }
# 
#    /** Solve A*X = B
#    @param  B   A Matrix with as many rows as A and any number of columns.
#    @return     X so that L*U*X = B(piv,:)
#    @exception  IllegalArgumentException Matrix row dimensions must agree.
#    @exception  RuntimeException  Matrix is singular.
#    */
# 
#    public Matrix solve (Matrix B) {
#       if (B.getRowDimension() != m) {
#          throw new IllegalArgumentException("Matrix row dimensions must agree.");
#       }
#       if (!this.isNonsingular()) {
#          throw new RuntimeException("Matrix is singular.");
#       }
# 
#       // Copy right hand side with pivoting
#       int nx = B.getColumnDimension();
#       Matrix Xmat = B.getMatrix(piv,0,nx-1);
#       double[][] X = Xmat.getArray();
# 
#       // Solve L*Y = B(piv,:)
#       for (int k = 0; k < n; k++) {
#          for (int i = k+1; i < n; i++) {
#             for (int j = 0; j < nx; j++) {
#                X[i][j] -= X[k][j]*LU[i][k];
#             }
#          }
#       }
#       // Solve U*X = Y;
#       for (int k = n-1; k >= 0; k--) {
#          for (int j = 0; j < nx; j++) {
#             X[k][j] /= LU[k][k];
#          }
#          for (int i = 0; i < k; i++) {
#             for (int j = 0; j < nx; j++) {
#                X[i][j] -= X[k][j]*LU[i][k];
#             }
#          }
#       }
#       return Xmat;
#    }
# end
# 
#-------------------------------------------------- 


#--------------------------------------------------
# #m = Matrix[[1,2,3],[4,5,6],[7,8,9]]
# m = Matrix[[1,0,0],[1,2,0],[1,4,8]]
# 
# lu = LUDecomposition.new(m)
# 
# puts lu.get_l
# puts
# puts lu.get_u
# puts
# puts lu.get_p
# puts
# 
# puts "p*l*u"
# puts lu.get_p * lu.get_l * lu.get_u
# puts
# puts "D = #{lu.det}"
#-------------------------------------------------- 
