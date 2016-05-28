module Algebra
class SVDecomposition 

	def signf(a,b)
		return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a)
	end

	def  pythag(a,b)
		absa = a.abs
		absb = b.abs
		if absa > absb
			c = absa * Math::sqrt(1.0+(absb/absa)**2)
		else
			c = (absb == 0.0 ? 0.0 : absb * Math::sqrt(1.0+(absa/absb)**2))
		end
		return c
	end

	def initialize(mat_a) # a : Matrix
		#: Calculates A = U * S * V.transpose

		# Working matrix, contains U at the end
		a = mat_a.deep_copy
		m = a.m
		n = a.n
		# W, as vector of diagonal emenents
		svec = Vector.of_size(n)
		# V, n by n matrix
		v = Matrix.zero(n,n)

		rv1 = Vector.of_size(n)
		g = 0.0
		scale = 0.0
		anorm = 0.0

		l = 0

		# Householder transformation to bidiagonal form
		n.times{|i|
			l = i + 2
			rv1[i]=scale*g
			g = s = scale = 0.0
			if i < m
				k = i
				while k < m do
					scale += a[k,i].abs
					k += 1
				end
				if scale != 0.0
					k = i
					while k < m do
						a[k,i] /= scale
						s += a[k,i] * a[k,i]
						k += 1
					end
					f=a[i,i]
					g = -signf(s**0.5,f)
					h=f*g-s
					a[i,i]=f-g
					j = l - 1
					while j < n do
						s = 0.0
						k = i
						while k < m do
							s += a[k,i]*a[k,j]
							k += 1
						end
						f=s/h
						k = i
						while k < m do
							a[k,j] += f*a[k,i]
							k += 1
						end
						j += 1
					end
					k = i
					while k < m do
						a[k,i] *= scale
						k += 1
					end
				end
			end
			svec[i] = scale * g
			g = s = scale = 0.0
			if i+1 <= m && i != n
				k=l-1
				while k < n do
					scale += a[i,k].abs
					k += 1
				end
				if scale != 0
					k= l -1
					while k < n do
						a[i,k] /= scale
						s += a[i,k]*a[i,k]
						k += 1
					end
					f=a[i,l-1]
					g = -signf(s**0.5,f)
					h=f*g-s
					a[i,l-1]=f-g
					k=l-1
					while k < n do
						rv1[k]=a[i,k]/h
						k += 1
					end
					j=l-1
					while j < m do
						s = 0.0
						k=l-1
						while k < n do
							s += a[j,k]*a[i,k]
							k += 1
						end
						k=l-1
						while k < n do
							a[j,k] += s*rv1[k]
							k += 1
						end
						j += 1
					end
					k=l-1
					while k < n do
						a[i,k] *= scale
						k += 1
					end
				end
			end
			anorm = [anorm, svec[i].abs + rv1[i].abs].max
		}

		# Accumulation of right-hand transformation 
		i = n - 1
		while i >= 0 do
			if i < n - 1
				if g != 0
					j = l
					while j < n do
						v[j,i]=(a[i,j]/a[i,l])/g
						j += 1
					end
					j = l
					while j < n do
						s = 0.0
						k = l
						while k < n do
							s += a[i,k]*v[k,j]
							k += 1
						end
						k = l
						while k < n do
							v[k,j] += s*v[k,i]
							k += 1
						end
						j += 1
					end
				end
				j = l
				while j < n do
					v[i,j] = v[j,i] = 0.0
					j += 1
				end
			end
			v[i,i]=1.0
			g=rv1[i]
			l=i
			i -= 1
		end

		# Accumulation of left-hand transformations
		i = [m,n].min - 1
		while i >= 0 do
			l=i+1
			g=svec[i]
			j=l
			while j<n do
				a[i,j]=0.0
				j += 1
			end
			if g != 0.0
				g=1.0/g
				j=l
				while j<n do
					s = 0.0
					k=l
					while k<m do
						s += a[k,i]*a[k,j]
						k += 1
					end
					f=(s/a[i,i])*g
					k=i
					while k<m do
						a[k,j] += f*a[k,i]
						k += 1
					end
					j += 1
				end
				j = i
				while j < m do
					a[j,i] *= g
					j += 1
				end
			else
				j = i
				while j < m do
					a[j,i]=0.0
					j += 1
				end
			end
			a[i,i] += 1.0
			i -= 1
		end

		# Diagonalization of the bidiagonal matrix
		k=n-1
		while k>=0 do
			# Do up to 30 iterations
			30.times {|its|
				flag=true
				l=k
				while l>=0 do
					nm=l-1
					if rv1[l].abs + anorm == anorm
						flag=false
						break
					end
					break if svec[nm].abs + anorm == anorm
					l -= 1
				end
				if flag
					c=0.0
					s=1.0
					i=l-1
					while i<k+1
						f=s*rv1[i]
						rv1[i]=c*rv1[i]
						break if f.abs + anorm == anorm
						g=svec[i]
						h=pythag(f,g)
						svec[i]=h
						h=1.0/h
						c=g*h
						s = -f*h
						m.times{|j|
							y=a[j,nm]
							z=a[j,i]
							a[j,nm]=y*c+z*s
							a[j,i]=z*c-y*s
						}
						i += 1
					end
				end
				z=svec[k]
				if l == k
					if z < 0.0
						svec[k] = -z
						n.times {|j| v[j,k] = -v[j,k]}
					end
					break
				end
				raise "Iterative SVD procedure not converged" if its == 29
				x=svec[l]
				nm=k-1
				y=svec[nm]
				g=rv1[nm]
				h=rv1[k]
				f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
				g = pythag(f,1.0)
				f=((x-z)*(x+z)+h*((y/(f+signf(g,f)))-h))/x
				c=s=1.0
				j=l
				while j<=nm do
					i=j+1
					g=rv1[i]
					y=svec[i]
					h=s*g
					g=c*g
					z = pythag(f,h)
					rv1[j]=z
					c=f/z
					s=h/z
					f=x*c+g*s
					g=g*c-x*s
					h=y*s
					y *= c
					n.times {|jj|
						x=v[jj,j]
						z=v[jj,i]
						v[jj,j]=x*c+z*s
						v[jj,i]=z*c-x*s
					}
					z = pythag(f,h)
					svec[j]=z
					if z != 0.0
						z=1.0/z
						c=f*z
						s=h*z
					end
					f=c*g+s*y
					x=c*y-s*g
					m.times {|jj|
						y=a[jj,j]
						z=a[jj,i]
						a[jj,j]=y*c+z*s
						a[jj,i]=z*c-y*s
					}
					j += 1
				end
				rv1[l]=0.0
				rv1[k]=f
				svec[k]=x
			}
			k -= 1
		end
		@u = a
		@v = v
		@s = svec

		@m = m
		@n = n
		return nil
	end

	def get_u
		#: Return the left singular vectors
		x = Matrix.zero(@m,@m)
		@m.times{|i|
			[@m,@n].min.times{|j|
				x[i,j] = @u[i,j]
			}
		}
		return x
	end

	def get_v
		#: Return the right singular vectors, transposed
		return @v.transpose
	end

	def get_s
		x = Matrix.zero(@m,@n)
		min = [@m,@n].min
		min.times{|i|
			x[i,i] = @s[i]
		}
		return x
	end

end
end

