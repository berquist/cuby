require "classes/objects/hessian.rb"

class Gradient < Array

	def Gradient.zero(natom)
		g = Gradient.new
		natom.times{g << Coordinate.new}
		return g
	end

	def Gradient.from_matrix(matrix)
		g = Gradient.new
		matrix.m.times{|i|
			g[i] = Coordinate[matrix[i,0], matrix[i,1], matrix[i,2]]
		}
		return g
	end

	def Gradient.from_vector(vector)
		g = Gradient.new
		raise "Wrong size of vector" unless vector.size % 3 == 0
		(vector.size/3).times{|i|
			g[i] = Coordinate[vector[i*3], vector[i*3+1], vector[i*3+2]]
		}
		return g
	end

	def Gradient.from_calculation_coordinates(geo_cuby, geo_new, grad_new)
		# Finds an orthogonal transformation matrix (not necessarily only a rotation) that maps
		# the calculation geometry onto the input one, and applies this transformation to
		# the gradient. The solution uses:
		# http://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem

		com_i = geo_cuby.center_of_mass
		com_o = geo_new.center_of_mass

		gi = geo_cuby.translate(-com_i)
		go = geo_new.translate(-com_o)

		p = gi.to_matrix
		q = go.to_matrix
		a = p.transpose * q
		u, s, v = a.svd

		# This will constrain the solution only to rotations, based on:
		# http://en.wikipedia.org/wiki/Kabsch_algorithm
		# (the form of the signmat matrix might be wrong here, see the link above)
		#--------------------------------------------------
		# sign =(v*u.transpose).det
		# sign = sign/sign.abs
		# signmat = Matrix.unit(3)
		# signmat[2,2] = sign
		#-------------------------------------------------- 


		rotm = v.transpose * u.transpose

		#Cuby.log.puts_debug "Coordinate transformation: largest error #{(p - q*rotm).max_abs} A"

		return Gradient.from_matrix(grad_new.to_matrix * rotm)
	end

	def to_s(format = "%15.8f")
		return self.map{|c| sprintf(format*3, c.x, c.y, c.z)}.join("\n")
	end

	def to_matrix
		m = Matrix.zero(self.size, 3)
		each_index{|i|
			crd = at(i)
			m[i,0], m[i,1], m[i,2] = crd.x, crd.y, crd.z
		}
		return m
	end

	def plus!(gradient)
		each_index{|i|
			at(i).plus!(gradient.at(i))
		}
	end

	def +(g2)
		raise(ArgumentError, "Gradients must be of same size") unless self.size == g2.size
		result = Gradient.new
		each_index{|i|
			result[i] = at(i) + g2.at(i)
		}
		return result
	end

	def -(g2)
		raise(ArgumentError, "Gradients must be of same size") unless self.size == g2.size
		result = Gradient.new
		each_index{|i|
			result[i] = at(i) - g2.at(i)
		}
		return result
	end

	def *(number)
		raise(ArgumentError, "Gradient can be multiplied only by a number") unless number.kind_of?(Numeric)
		result = Gradient.new
		each_index{|i|
			result[i] = at(i) * number
		}
		return result
	end

	def append(g2)
		g = Gradient.new
		each{|c| g << c}
		g2.each{|c| g << c}
		return g
	end

	def project_out_trans_rot!(geometry)
		projector = Hessian.transrot_removal_projector(geometry)
		v = projector * self.to_vector
		self.update_from_vector!(v)
		return nil
	end

	#=======================================================================
	# Vector access
	#=======================================================================
	
	def to_vector # => Vector
		array = []
		each {|crd|
			array << crd.x
			array << crd.y
			array << crd.z
		}
		return Vector.from_array(array)
	end

	def update_from_vector!(vector) # => nil
		# test size
		raise "Wrong size of vector" unless vector.size == size * 3
		# copy values
		each_index {|i|
			crd = at(i)
			crd.x = vector[i*3]
			crd.y = vector[i*3+1]
			crd.z = vector[i*3+2]
		}
		return nil
	end

	#=======================================================================
	# Misc
	#=======================================================================
	
	def rms
		sumsq = 0.0
		each{|c| sumsq += c.x**2 + c.y**2 + c.z**2 }
		sumsq /= size * 3
		return sumsq**0.5
	end

	def abs
		sumsq = 0.0
		each{|c| sumsq += c.x**2 + c.y**2 + c.z**2 }
		return sumsq**0.5
	end

	def max_component
		max = 0.0
		each{|c|
			c.each{|v|
				max = v.abs if v.abs > max
			}
		}
		return max
	end

	#=======================================================================
	# Gradient of selection
	#=======================================================================
	
	def for_atomlist(atomlist)
		g = Gradient.new
		atomlist.each{|i|
			g << self.at(i)
		}
		return g
	end

	#=======================================================================
	# Gradient from internal coordinates
	#=======================================================================
	
	def add_internal_dist(geometry, grad, at_i, at_j)
		u = geometry.at(at_i) - geometry.at(at_j) # Atom-atom vector
		gcart = grad / u.abs * u
		self.at(at_i).plus!(gcart)
		self.at(at_j).minus!(gcart)
		return nil
	end

	def add_internal_angle(geometry, grad, at_i, at_j, at_k)
		u = geometry.at(at_i) - geometry.at(at_j)
		v = geometry.at(at_k) - geometry.at(at_j)
		ru = u.abs
		rv = v.abs
		dot = u.dot(v)

		g_i = grad / (1.0 - dot**2 / ru**2 / rv**2)**0.5 * -(v/ru/rv - u*dot/ru**3/rv)
		g_k = grad / (1.0 - dot**2 / ru**2 / rv**2)**0.5 * -(u/ru/rv - v*dot/ru/rv**3)
		self.at(at_i).plus!(g_i)
		self.at(at_j).plus!(-g_i - g_k)
		self.at(at_k).plus!(g_k)
		return nil
	end

	def add_internal_torsion(geometry, grad, at_i, at_j, at_k, at_l)
		# Atoms:
		i, j, k, l = geometry.at(at_i), geometry.at(at_j), geometry.at(at_k), geometry.at(at_l)
		# Important vectors:
		u = j - i
		v = k - j
		w = l - k
		t = i - k
		s = j - l
		v_abs = v.abs

		# No gradient if one of the distances is 0
		[u,v,w,t,s].each{|vec|
			return nil if vec.abs == 0.0
		}

		# Linear case check
		a1 = u.angle(v)
		a2 = v.angle(w)
		if a1 == 0.0 || a1 == Math::PI ||
		   a2 == 0.0 || a2 == Math::PI
			return nil
		end

		cuv = u.cross_product(v)
		cvw = v.cross_product(w)

		ex1 = v_abs * u.dot(cvw) 
		ex2 = cuv.dot(cvw) 

		denom = (ex1**2+ex2**2)

		self.at(at_i).plus!( 
			(v.cross_product(cvw)*ex1 - cvw*v_abs*ex2) * grad / denom
		)
		self.at(at_l).plus!( 
			(v.cross_product(cuv)*ex1 + cuv*v_abs*ex2) * grad / denom
		)
		self.at(at_j).plus!( 
			(-u.dot(cvw)/v_abs*ex2*v + cvw*v_abs*ex2 + u.cross_product(w)*v_abs*ex2 +
			t.cross_product(cvw)*ex1 - cuv.cross_product(w)*ex1) * grad / denom
		)
		self.at(at_k).plus!( 
			(u.dot(cvw)/v_abs*ex2*v + u.cross_product(s)*v_abs*ex2 +
			u.cross_product(cvw)*ex1 - cuv.cross_product(s)*ex1) * grad / denom
		)
		return nil
	end
end
