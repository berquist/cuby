module Algebra
class QRDecomposition

	def initialize(a)
		@qr = a.to_a
		@m = a.m
		@n = a.n

		@rdiag = Array.new(@n)

		@n.times {|k|
			# Compute 2-norm of k-th column without under/overflow
			nrm = 0
			(@m-k).times{|ii| # for (int i = k; i < m; i++)
				i = ii + k
				nrm = (nrm**2 + @qr[i][k]**2)**0.5;
			}

			if nrm != 0.0
				# Form k-th Householder vector
				if @qr[k][k] < 0
					nrm = -nrm
				end
				(@m-k).times{|ii| # (int i = k; i < m; i++)
					i = ii + k
					@qr[i][k] /= nrm
				}
				@qr[k][k] += 1.0

				# Apply transformation to remaining columns
				(@n-k-1).times { |jj| # for (int j = k+1; j < n; j++)
					j = jj + k + 1
					s = 0.0
					(@m-k).times{|ii| # (int i = k; i < m; i++)
						i = ii + k
						s += @qr[i][k]*@qr[i][j]
					}
					s = -s/@qr[k][k]
					(@m-k).times{|ii| # (int i = k; i < m; i++)
						i = ii + k
						@qr[i][j] += s*@qr[i][k]
					}
				}
			end
			@rdiag[k] = -nrm
		}
	end

	def is_full_rank?
		# true if R, and hence A, has full rank
		@n.times {|j|
			return false if @rdiag[j] == 0
		}
		return true
	end

	def get_h # => Matrix
		# Return the Householder vectors
		# Lower trapezoidal matrix whose columns define the reflections
		h = Matrix.zero(@m,@n)
		@m.times {|i|
			@n.times {|j|
				if i >= j
					h[i,j] = @qr[i][j]
				#else
				#	h[i,j] = 0.0
				end
			}
		}
		return h
	end

	def get_r # => Matrix
		# Return the upper triangular factor
		r = Matrix.zero(@m,@n)
		@m.times {|i|
			@n.times {|j|
				if i < j
					r[i,j] = @qr[i][j]
				elsif i == j
					r[i,j] = @rdiag[i]
				else
					#r[i,j] = 0.0
				end
			}
		}
		return r
	end

	def get_q
		# Generate and return the (economy-sized) orthogonal factor
		q = Matrix.zero(@m,@m)
		minmn = [@m,@n].min
		minmn.times{|kk| # for (int k = n-1; k >= 0; k--)
			k = minmn - kk -1
			@m.times {|i|
				q[i,k] = 0.0
			}
			q[k,k] = 1.0
			(minmn-k).times {|jj| #  for (int j = k; j < n; j++)
				j = jj + k
				if @qr[k][k] != 0
					s = 0.0
					(@m-k).times {|ii| # for (int i = k; i < m; i++)
						i = ii + k
						s += @qr[i][k]*q[i,j]
					}
					s = -s/@qr[k][k]
					(@m-k).times {|ii| # for (int i = k; i < m; i++)
						i = ii + k
						q[i,j] += s*@qr[i][k]
					}
				end
			}
		}
		return q
	end

	def solve(b)
		# b i s Matrix, with as many rows as A and any number of columns
		# returns X that minimizes the two norm of Q*R*X-B

		raise "Matrix row dimensions must agree." if b.m != @m
		raise "Matrix is rank deficient" unless is_full_rank?

		# Copy right hand side
		nx = b.n
		x = b.dup

		# Compute Y = transpose(Q)*B
		@n.times {|k|
			nx.times {|j|
				s = 0.0
				(@m-k).times {|ii| # for (int i = k; i < m; i++)
					i = ii + k
					s += @qr[i][k]*x[i,j]
				}
				s = -s/@qr[k][k]
				(@m-k).times {|ii| # for (int i = k; i < m; i++)
					i = ii + k
					x[i,j] += s*@qr[i][k]
				}
			}
		}

		# Solve R*X = Y
		@n.times {|kk| # for (int k = n-1; k >= 0; k--)
			k = @n - k - 1
			nx.times {|j|
				x[k,j] /= @rdiag[k]
			}
			k.times {|i|
				nx.times {|j|
					x[i,j] -= x[k,j]*@qr[i][k]
				}
			}
		}
		#return x.submatrix(0,n-1,0,nx-1) originaly, by indexes
		return x.submatrix(0,0,n,nx) #new, by size
	end
end
end

