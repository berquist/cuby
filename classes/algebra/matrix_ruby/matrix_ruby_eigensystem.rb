module Algebra
class EigenvalueDecomposition
	attr_reader :issymmetric
	# variables @n, @issymmetric, @d, @e, @v, @h, @ort
	
	def initialize(matrix)
		a = matrix.to_a
		@n = matrix.n
		@v = [] # will be n x n
		@d = Array.new(@n,0.0)
		@e = Array.new(@n,0.0)

		@issymmetric = true
		j = 0; i = 0
		while (j < @n) && @issymmetric
			while (i < @n) && @issymmetric
				 @issymmetric = (a[i][j] == a[j][i])
				 i += 1
			end
			j += 1
		end

		if @issymmetric
			@n.times{|i|
				@v[i] = []
				@n.times{|j|
					 @v[i][j] = a[i][j]
				}
			}

			# Tridiagonalize
			tred2

			# Diagonalize
			tql2
		else
			@h = []
			@ort = Array.new(@n,0.0)
			@n.times{|i|
				@h[i] = []
				@v[i] = [] ###
				@n.times{|j|
					 @h[i][j] = a[i][j]
					 @v[i][j] = 0.0 ### build the array
				}
			}

			# Reduce to Hessenberg form
			orthes

			#Reduce Hessenberg to real Schur form
			hqr2
		end
	end

	def get_d
		x = Matrix.zero(@n,@n)
		@n.times{|i|
			x[i,i] = @d[i]
			if (@e[i] > 0)
				x[i,i+1] = @e[i]
			elsif (@e[i] < 0)
				x[i,i-1] = @e[i]
			end
		}
		return x
	end

	def get_v
		return Matrix.from_array(@v)
	end

	def get_real_eigval
		return @d
	end

	def get_imag_eigval
		return @e
	end

#    /** Return the real parts of the eigenvalues
#    @return     real(diag(D))
#    */
# 
#    public double[] getRealEigenvalues () {
#       return d;
#    }
# 
#    /** Return the imaginary parts of the eigenvalues
#    @return     imag(diag(D))
#    */
# 
#    public double[] getImagEigenvalues () {
#       return e;
#    }
	
	def tred2
		# Symmetric Householder reduction to tridiagonal form

		@n.times{|j|
			@d[j] = @v[@n-1][j]
		}

		@n.times{|ii|
			i = @n - 1 - ii # @n-1 downto 0

			# Scale to avoid under/overflow
			scale = 0.0
			h = 0.0
			i.times{|k|
				scale += @d[k].abs
			}
			if scale == 0.0
				@e[i] = @d[i-1]
				i.times{|j|
					@d[j] = @v[i-1][j]
					@v[i][j] = 0.0
					@v[j][i] = 0.0
				}
			else
				# Generate Householder vector
				i.times{|k|
					@d[k] /= scale
					h += @d[k] * @d[k]
				}
				f = @d[i-1]
				g = Math.sqrt(h)
				g = -g if (f > 0.0)
				@e[i] = scale * g
				h = h - f * g
				@d[i-1] = f - g
				i.times{|j|
					@e[j] = 0.0
				}
				# Apply similarity transformation to remaining columns
				i.times{|j|
					f = @d[j]
					@v[j][i] = f
					g = @e[j] + @v[j][j] * f
					(i-j-1).times {|kk| # for (int k = j+1; k <= i-1; k++)
						k = j+1+kk
						g += @v[k][j] * @d[k]
						@e[k] += @v[k][j] * f
					}
					@e[j] = g
				}
				f = 0.0
				i.times{|j|
					@e[j] /= h
					f += @e[j] * @d[j]
				}
				hh = f / (h + h)
				i.times{|j|
					@e[j] -= hh * @d[j]
				}
				i.times{|j|
					f = @d[j]
					g = @e[j]
					(i-j).times {|kk| # for (int k = j; k <= i-1; k++)
						 k = j+kk
						 @v[k][j] -= (f * @e[k] + g * @d[k])
					}
					@d[j] = @v[i-1][j]
					@v[i][j] = 0.0
				}
			end
			@d[i] = h
		}

		# Accumulate transformations
		(@n-1).times {|i|
			@v[@n-1][i] = @v[i][i]
			@v[i][i] = 1.0
			h = @d[i+1] 
			if (h != 0.0)
				(i+1).times{|k| # for (int k = 0; k <= i; k++)
					@d[k] = @v[k][i+1] / h
				}
				(i+1).times{|j| # for (int j = 0; j <= i; j++)
					g = 0.0
					(i+1).times{|k| # for (int k = 0; k <= i; k++)
						g += @v[k][i+1] * @v[k][j]
					}
					(i+1).times{|k| # for (int k = 0; k <= i; k++)
						@v[k][j] -= g * @d[k]
					}
				}
			end
			(i+1).times{|k| # for (int k = 0; k <= i; k++)
				 @v[k][i+1] = 0.0;
			}
		}
		@n.times{|j|
			@d[j] = @v[@n-1][j]
			@v[@n-1][j] = 0.0
		}
		@v[@n-1][@n-1] = 1.0
		@e[0] = 0.0
	end

	def tql2
		(@n-1).times {|ii|
			i = ii+1
			@e[i-1] = @e[i]
		}
		@e[@n-1] = 0.0

		f = 0.0
		tst1 = 0.0
		eps = 2.0**-52
		
		@n.times {|l|
			# Find small subdiagonal element
			tst1 = [tst1, @d[l].abs + @e[l].abs].max
			m = l
			while (m < @n)
				break if @e[m].abs <= eps*tst1
				m += 1
			end

			# If m == l, @d[l] is an eigenvalue, otherwise, iterate
			if (m > l)
				iter = 0
				do_cycle = true
				while do_cycle
					iter += 1
					# Compute implicit shift
					g = @d[l]
					p = (@d[l+1] - g) / (2.0 * @e[l])
					r = (p**2+1.0)**0.5 # euclidean distance
					r = -r if p < 0.0
					@d[l] = @e[l] / (p + r)
					@d[l+1] = @e[l] * (p + r)
					dl1 = @d[l+1]
					h = g - @d[l] 
					(@n-l-2).times {|ii| # for (int i = l+2; i < @n; i++)
						i = ii + l + 2
						@d[i] -= h
					}
					f = f + h

					# Implicit QL transformation
					p = @d[m]
					c = 1.0
					c2 = c
					c3 = c
					el1 = @e[l+1]
					s = 0.0
					s2 = 0.0
					(m-l).times{ |ii|# for (int i = m-1; i >= l; i--) # tested
						i = m-1-ii
						c3 = c2
						c2 = c
						s2 = s
						g = c * @e[i]
						h = c * p
						r = (p**2 + @e[i]**2)**0.5
						@e[i+1] = s * r
						s = @e[i] / r
						c = p / r
						p = c * @d[i] - s * g
						@d[i+1] = h + s * (c * g + s * @d[i])
						# Accumulate transformation
						@n.times{|k|
							h = @v[k][i+1]
							@v[k][i+1] = s * @v[k][i] + c * h
							@v[k][i] = c * @v[k][i] - s * h
						}
					}
					p = -s * s2 * c3 * el1 * @e[l] / dl1
					@e[l] = s * p
					@d[l] = c * p

					# check convergence
					do_cycle = @e[l].abs > eps*tst1
				end # while
			end
			@d[l] = @d[l] + f
			@e[l] = 0.0
		}

		# sort eigenvalues and corresponding vectors
		(@n-1).times{|i|
			k = i
			p = @d[i]
			(@n-i-1).times {|jj| #for (int j = i+1; j < @n; j++)
				j = jj + i + 1
				if (@d[j] < p)
					k = j
					p = @d[j]
				end
			}
			if (k != i)
				@d[k] = @d[i]
				@d[i] = p
				@n.times{|j|
					p = @v[j][i]
					@v[j][i] = @v[j][k]
					@v[j][k] = p
				}
			end
		}
	end

	# Nonsymmetric reduction to Hessenberg form
	def orthes
		low = 0
		high = @n-1
		(high - low - 1).times {|mm| #for (int m = low+1; m <= high-1; m++)
			m = mm + low + 1
			#  Scale column
			scale = 0.0
			(high - m).times {|ii| # for (int i = m; i <= high; i++)
				i = ii + m
				scale = scale + @h[i][m-1].abs
			}
			if scale != 0.0
				# Compute Householder transformation
				h = 0.0
				(high - m + 1).times {|ii| # for (int i = high; i >= m; i--) # tested
					i = high - ii
					@ort[i] = @h[i][m-1]/scale
					h += @ort[i] * @ort[i]
				}
				g = Math.sqrt(h)
				g = -g if (@ort[m] > 0)
				h = h - @ort[m] * g
				@ort[m] = @ort[m] - g

				# Apply Householder similarity transformation
				(@n-m).times {|jj| # for (int j = m; j < @n; j++) # tested
					j = jj + m
					f = 0.0
					(high - m + 1).times {|ii| # for (int i = high; i >= m; i--) # tested
						i = high - ii
						f += @ort[i]*@h[i][j]
					}
					f = f/h
					(high - m + 1).times {|ii| # for (int i = m; i <= high; i++)
						i = ii + m
						@h[i][j] -= f*@ort[i]
					}
				}

				(high + 1).times {|i|
					f = 0.0
					(high - m + 1).times {|jj| # for (int i = high; i >= m; i--) # tested
						j = high - jj
						f += @ort[j]*@h[i][j]
					}
					f = f/h
					(high - m + 1).times {|jj| # for (int i = m; i <= high; i++)
						j = jj + m
						@h[i][j] -= f*@ort[j]
					}
				}
				@ort[m] = scale*@ort[m]
				@h[m][m-1] = scale*g
			end
		}

		# Accumulate transformations
		@n.times{|i|
			@n.times{|j|
				 @v[i][j] = (i == j ? 1.0 : 0.0) #???
			}
		}
		(high-1-low).times {|mm| # for (int m = high-1; m >= low+1; m--) # tested
			m = high - 1 - mm
			if (@h[m][m-1] != 0.0)
				(high - m).times {|ii| # for (int i = m+1; i <= high; i++) # tested
					i = ii + m + 1
					@ort[i] = @h[i][m-1]
				}
				(high + 1 - m).times {|jj| # for (int j = m; j <= high; j++) # tested
					j = jj + m
					g = 0.0
					(high + 1 - m).times {|ii|
						i = ii + m
						g += @ort[i] * @v[i][j]
					}
					# Double division avoids possible underflow
					g = (g / @ort[m]) / @h[m][m-1]
					(high + 1 - m).times {|ii|
						i = ii + m
						@v[i][j] += g * @ort[i]
					}
				}
			end
		}
	end

	# returns [cdivr,cdivi;
	def  cdiv(xr, xi, yr, yi)
		if (yr.abs > yi.abs)
			r = yi/yr
			d = yr + r*yi
			return [(xr + r*xi)/d, (xi - r*xr)/d]
		else
			r = yr/yi
			d = yi + r*yr
			return [(r*xr + xi)/d, (r*xi - xr)/d]
		end
	end
	private :cdiv

	# Nonsymmetric reduction from Hessenberg to real Schur form
	def hqr2
		# initialize
		nn = @n
		@n = nn-1
		low = 0
		high = nn-1
		eps = 2.0**-52
		exshift = 0.0
		p=0.0; q=0.0; r=0.0; s=0.0; z=0.0
		t=0.0; w=0.0; x=0.0; y=0.0

		# Store roots isolated by balanc and compute matrix norm
		#-------------------------------------------------- 
		norm = 0.0
		nn.times{|i|
			if (i < low || i > high)
				@d[i] = @h[i][i]
				@e[i] = 0.0
			end
			jstart = [i-1,0].max
			(nn-jstart).times {|jj| #  for (int j = MathExt.max(i-1,0); j < nn; j++) # tested
				j = jstart + jj
				norm = norm + @h[i][j].abs
			}
		}

		# Outer loop over eigenvalue index
		iter = 0
		while (@n >= low)
			# Look for single small sub-diagonal element
			l = @n
			while (l > low)
				s = @h[l-1][l-1].abs + @h[l][l].abs
				if (s == 0.0)
					s = norm
				end
				if (@h[l][l-1].abs < eps * s)
					break
				end
				l -= 1
			end
			# Check for convergence
			if (l == @n)
				# One root found
				@h[@n][@n] = @h[@n][@n] + exshift
				@d[@n] = @h[@n][@n]
				@e[@n] = 0.0
				@n -= 1
				iter = 0
			elsif (l == @n-1)
				# Two roots found
				w = @h[@n][@n-1] * @h[@n-1][@n]
				p = (@h[@n-1][@n-1] - @h[@n][@n]) / 2.0
				q = p * p + w
				z = Math.sqrt(q.abs)
				@h[@n][@n] = @h[@n][@n] + exshift
				@h[@n-1][@n-1] = @h[@n-1][@n-1] + exshift
				x = @h[@n][@n]

				if (q >= 0) # real pair
					if (p >= 0)
						z = p + z
					else
						z = p - z
					end
					@d[@n-1] = x + z
					@d[@n] = @d[@n-1]
					if (z != 0.0)
						@d[@n] = x - w / z
					end
					@e[@n-1] = 0.0
					@e[@n] = 0.0
					x = @h[@n][@n-1]
					s = x.abs + z.abs
					p = x / s
					q = z / s
					r = Math.sqrt(p * p+q * q)
					p = p / r
					q = q / r

					# row modification
					(nn-@n+1).times {|jj| # for (int j = @n-1; j < nn; j++)
						j = jj +  @n-1
						z = @h[@n-1][j]
						@h[@n-1][j] = q * z + p * @h[@n][j]
						@h[@n][j] = q * @h[@n][j] - p * z
					}

					# Column modification
					(@n+1).times {|i| #  for (int i = 0; i <= @n; i++)
						z = @h[i][@n-1]
						@h[i][@n-1] = q * z + p * @h[i][@n]
						@h[i][@n] = q * @h[i][@n] - p * z
					}

					# Accumulate transformations
					(high-low+1).times {|ii| # for (int i = low; i <= high; i++)
					       i = ii + low
					       z = @v[i][@n-1]
					       @v[i][@n-1] = q * z + p * @v[i][@n]
					       @v[i][@n] = q * @v[i][@n] - p * z
					}
				else # Complex pair
					@d[@n-1] = x + p
					@d[@n] = x + p
					@e[@n-1] = z
					@e[@n] = -z
				end
				@n = @n - 2
				iter = 0
			else # No convergence yet
				# Form shift
				x = @h[@n][@n]
				y = 0.0
				w = 0.0
				if (l < @n)
					y = @h[@n-1][@n-1]
					@h[@n][@n-1] * @h[@n-1][@n]
				end

				# Wilkinson''s original ad hoc shift
				if (iter == 10)
					exshift += x
					(@n+1-low).times {|ii| # for (int i = low; i <= @n; i++)
						i = ii + low
						 @h[i][i] -= x
					}
					s = @h[@n][@n-1].abs + @h[@n-1][@n-2].abs
					x = y = 0.75 * s
					w = -0.4375 * s * s
				end

				# MATLAB''s new ad hoc shift
				if (iter == 30)
					s = (y - x) / 2.0
					s = s * s + w
					if (s > 0)
						s = Math.sqrt(s)
						if (y < x)
							s = -s
						end
						s = x - w / ((y - x) / 2.0 + s)
						(@n+1-low).times {|ii| # for (int i = low; i <= @n; i++)
							i = ii + low
							@h[i][i] -= s
						}
						exshift += s
						x = y = w = 0.964
					end
				end

				iter = iter + 1 # Could check iteration count here

				# Look for two consecutive small sub-diagonal elements
				m = @n-2
				while (m >= l)
					z = @h[m][m]
					r = x - z
					s = y - z
					p = (r * s - w) / @h[m+1][m] + @h[m][m+1]
					q = @h[m+1][m+1] - z - r - s
					r = @h[m+2][m+1]
					s = p.abs + q.abs + r.abs
					p = p / s
					q = q / s
					r = r / s
					break if (m == l)
					break if (@h[m][m-1].abs * (q.abs + r.abs) < eps * (p.abs * (@h[m-1][m-1].abs + z.abs + @h[m+1][m+1].abs)))
					m -= 1
				end

				(@n+1-m-2).times {|ii| #for (int i = m+2; i <= @n; i++)
					i = m+2+ii
					@h[i][i-2] = 0.0
					@h[i][i-3] = 0.0 if (i > m+2)
				}

				# Double QR step involving rows l:@n and columns m:n
				(@n-m).times {|kk| # for (int k = m; k <= @n-1; k++)
					k = kk + m
					notlast = (k != @n-1)
					if (k != m)
						p = @h[k][k-1]
						q = @h[k+1][k-1]
						r = (notlast ? @h[k+2][k-1] : 0.0)
						x = p.abs + q.abs + r.abs
						if (x != 0.0)
							p = p / x
							q = q / x
							r = r / x
						end
					end
					break if (x == 0.0)
					s = Math.sqrt(p * p + q * q + r * r)
					s = -s if (p < 0)
					if (s != 0)
						if (k != m)
							@h[k][k-1] = -s * x
						elsif (l != m)
							@h[k][k-1] = -@h[k][k-1]
						end
						p = p + s
						x = p / s
						y = q / s
						z = r / s
						q = q / p
						r = r / p

						# Row modification
						(nn - k).times { |jj| # for (int j = k; j < nn; j++)
							j = jj + k
							p = @h[k][j] + q * @h[k+1][j]
							if notlast
								p = p + r * @h[k+2][j]
								@h[k+2][j] = @h[k+2][j] - p * z
							end
							@h[k][j] = @h[k][j] - p * x
							@h[k+1][j] = @h[k+1][j] - p * y
						}

						# Column modification
						([@n,k+3].min + 1).times {|i| # for (int i = 0; i <= MathExt.min(@n,k+3); i++)
							p = x * @h[i][k] + y * @h[i][k+1]
							if notlast
								p = p + z * @h[i][k+2]
								@h[i][k+2] = @h[i][k+2] - p * r
							end
							@h[i][k] = @h[i][k] - p
							@h[i][k+1] = @h[i][k+1] - p * q
						}

						# Accumulate transformations
						(high+1-low).times {|ii| # for (int i = low; i <= high; i++)
							i = ii+low
							p = x * @v[i][k] + y * @v[i][k+1]
							if (notlast)
								p = p + z * @v[i][k+2]
								@v[i][k+2] = @v[i][k+2] - p * r
							end
							@v[i][k] = @v[i][k] - p
							@v[i][k+1] = @v[i][k+1] - p * q
						}
					end
				}
			end
		end

		# Backsubstitute to find vectors of upper triangular form
		return if (norm == 0.0)

		nn.times {|nnc| # for (@n = nn-1; @n >= 0; @n--)
			@n = nn-1 -nnc
			p = @d[@n]
			q = @e[@n]
			# Real vector
			if (q == 0)
				l = @n
				@h[@n][@n] = 1.0
				@n.times {|ii| # for (int i = @n-1; i >= 0; i--)
					i = @n-1-ii
					w = @h[i][i] - p
					r = 0.0
					(@n+1-l).times {|jj| # for (int j = l; j <= @n; j++)
						j = jj+l
						r = r + @h[i][j] * @h[j][@n]
					}
					if (@e[i] < 0.0)
						z = w
						s = r
					else
						l = i
						if (@e[i] == 0.0)
							if (w != 0.0)
								@h[i][@n] = -r / w
							else
								@h[i][@n] = -r / (eps * norm)
							end
						else # Solve real equations
							x = @h[i][i+1]
							y = @h[i+1][i]
							q = (@d[i] - p) * (@d[i] - p) + @e[i] * @e[i]
							t = (x * s - z * r) / q
							@h[i][@n] = t
							if x.abs > z.abs
								@h[i+1][@n] = (-r - w * t) / x
							else
								@h[i+1][@n] = (-s - y * t) / z
							end
						end
						# Overflow control
						t = @h[i][@n].abs
						if ((eps * t) * t > 1)
							(@n+1-i).times {|jj| # for (int j = i; j <= @n; j++)
								j = jj + i
								@h[j][@n] = @h[j][@n] / t
							}
						end
					end
				}
			elsif (q < 0) # Complex vector
				l = @n-1
				if @h[@n][@n-1].abs > @h[@n-1][@n].abs
					@h[@n-1][@n-1] = q / @h[@n][@n-1]
					@h[@n-1][@n] = -(@h[@n][@n] - p) / @h[@n][@n-1]
				else
					@h[@n-1][@n-1], @h[@n-1][@n] = cdiv(0.0,-@h[@n-1][@n],@h[@n-1][@n-1]-p,q)
					#@h[@n-1][@n-1] = cdivr
					#@h[@n-1][@n] = cdivi
				end
				@h[@n][@n-1] = 0.0
				@h[@n][@n] = 1.0
				(@n-2+1).times {|ii| # for (int i = @n-2; i >= 0; i--)
					i = @n-2 - ii
					ra = 0.0
					sa = 0.0
					(@n+1-l).times {|jj| # for (int j = l; j <= @n; j++)
						j = jj+l
						ra = ra + @h[i][j] * @h[j][@n-1]
						sa = sa + @h[i][j] * @h[j][@n]
					}
					w = @h[i][i] - p
					if (@e[i] < 0.0)
						z = w
						r = ra
						s = sa
					else
						l = i
						if (@e[i] == 0)
							@h[i][@n-1], @h[i][@n] = cdiv(-ra,-sa,w,q)
							#@h[i][@n-1] = cdivr
							#@h[i][@n] = cdivi
						else
							# Solve complex equations
							x = @h[i][i+1]
							y = @h[i+1][i]
							vr = (@d[i] - p) * (@d[i] - p) + @e[i] * @e[i] - q * q
							vi = (@d[i] - p) * 2.0 * q
							if (vr == 0.0 && vi == 0.0)
								vr = eps * norm * (w.abs + q.abs + x.abs + y.abs + z.abs)
							end
							@h[i][@n-1], @h[i][@n] = cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi)
							#@h[i][@n-1] = cdivr
							#@h[i][@n] = cdivi
							if x.abs > z.abs + q.abs
								@h[i+1][@n-1] = (-ra - w * @h[i][@n-1] + q * @h[i][@n]) / x
								@h[i+1][@n] = (-sa - w * @h[i][@n] - q * @h[i][@n-1]) / x
							else
								@h[i+1][@n-1], @h[i+1][@n] = cdiv(-r-y*@h[i][@n-1],-s-y*@h[i][@n],z,q)
								#@h[i+1][@n-1] = cdivr
								#@h[i+1][@n] = cdivi
							end
						end

						# Overflow control
						t = [@h[i][@n-1].abs, @h[i][@n].abs].max
						if ((eps * t) * t > 1)
							(@n+1-i).times {|jj| # for (int j = i; j <= @n; j++)
								j = jj + i
								@h[j][@n-1] = @h[j][@n-1] / t
								@h[j][@n] = @h[j][@n] / t
							}
						end
					end
				}
			end
		}

		# Vectors of isolated roots
		nn.times {|i|
			if (i < low || i > high)
				(nn-i).times {|jj| # for (int j = i; j < nn; j++)
					j = jj + i
					@v[i][j] = @h[i][j]
				}
			end
		}

		# Back transformation to get eigenvectors of original matrix
		(nn-low).times {|jj| # for (int j = nn-1; j >= low; j--) tested
			j = nn - jj - 1
			(high+1 - low).times {|ii| # for (int i = low; i <= high; i++) # tested
				i = ii + low
				z = 0.0
				([j,high].min + 1 - low).times {|kk| # for (int k = low; k <= MathExt.min(j,high); k++)
					k = kk + low
					z = z + @v[i][k] * @h[k][j]
				}
				@v[i][j] = z
			}
		}
	end
end
end



################################################################################

#--------------------------------------------------
# m = Matrix[[1,2,3],[2,3,4],[3,4,5]]
# e = EigenvalueDecomposition.new(m)
# 
# puts "D"
# puts e.get_d
# puts
# puts "V"
# puts e.get_v
# puts
# puts "V*D*V'"
# puts e.get_v * e.get_d * e.get_v.transpose
# puts
# 
# puts "And now asymmetric one:"
# #m = Matrix[[1,2,3],[4,5,6],[7,8,9]]
# m = Matrix[[3,-2],[4,-1]]
# e = EigenvalueDecomposition.new(m)
# puts "D"
# puts e.get_d
# puts
# puts "V"
# puts e.get_v
# puts
# puts "real eigenvalues"
# puts e.get_real_eigval
# puts "imag eigenvalues"
# puts e.get_imag_eigval
#-------------------------------------------------- 
