################################################################################
#
# Module Jacobian
#
# Author: Jan Rezac
# Date created: 2009-05-08
# License: Cuby license
# Description: Jacobian matrix for internal coordinates
# Extends: N/A
# Status: Works, possible improvements
#
################################################################################

#: More detailed description

module Jacobian

	# Displacement for numerical derivatives
	@@num_dx = 0.0001 # Angstrom

	def Jacobian.signf(a,m,n)
		return  1.0 if a == m
		return -1.0 if a == n
		return  0.0
	end

	#=======================================================================
	# Coordinate derivatives
	#=======================================================================
	
	def Jacobian.dq_dx(type, cartesians, coord_i)
		case type
		when :distance
			return Jacobian.dq_dx_bond(cartesians, coord_i)
		when :distance_inv
			return Jacobian.dq_dx_distinv(cartesians, coord_i)
		when :angle
			return Jacobian.dq_dx_angle(cartesians, coord_i)
		when :torsion
			return Jacobian.dq_dx_torsion(cartesians, coord_i)
		else
			raise "Wrong internal coordinate type #{type}"
		end
	end

	def Jacobian.dq2_dxdy(type, cartesians, coord_i, coord_j)
		case type
		when :distance
			return Jacobian.dq2_dxdy_bond_num(cartesians, coord_i, coord_j)
		when :angle
			return Jacobian.dq2_dxdy_angle_num(cartesians, coord_i, coord_j)
		when :torsion
			return Jacobian.dq2_dxdy_torsion_num(cartesians, coord_i, coord_j)
		else
			raise "Wrong internal coordinate type #{type}"
		end
	end

	#=======================================================================
	# Elementary coordinates - first derivatives
	#=======================================================================

	def Jacobian.dq_dx_bond(cartesians, coord_i)
		#: Arguments:
		#* cartesians - Array (2) of Coordinate
		#* coord_i - index of the coordinate (0,1)

		u = (cartesians[0] - cartesians[1]).normalize
		sign = Jacobian.signf(coord_i, 0, 1)

		return sign * u
	end

	def Jacobian.dq_dx_distinv(cartesians, coord_i)
		#: Arguments:
		#* cartesians - Array (2) of Coordinate
		#* coord_i - index of the coordinate (0,1)
		
		u = -1.0 / (2.0 * (cartesians[0] - cartesians[1]).abs**2) * (cartesians[0] - cartesians[1]).normalize
		sign = Jacobian.signf(coord_i, 0, 1)
		
		return sign * u
	end

	def Jacobian.dq_dx_angle(cartesians, coord_i)
		#: Arguments:
		#* cartesians - Array (3) of Coordinate
		#* coord_i - index of the coordinate (0,1,2)

		u = cartesians[0] - cartesians[1]
		v = cartesians[2] - cartesians[1]

		lu = u.abs
		lv = v.abs

		u.normalize!
		v.normalize!

		# Get perpendicular vector	
		w = u.cross_product(v)
		if w == Coordinate[0.0,0.0,0.0] # u and v are parallel
			w = u.cross_product(Coordinate[1.0,-1.0,1.0])
		end
		if w == Coordinate[0.0,0.0,0.0] # u and the first guess are parallel, try another guess
			w = u.cross_product(Coordinate[-1.0,1.0,1.0])
		end
		w.normalize!

		# Signs
		sign1 = Jacobian.signf(coord_i, 0, 1)
		sign2 = Jacobian.signf(coord_i, 2, 1)

		# Derivative
		return 	sign1 / lu * u.cross_product(w) + 
			sign2 / lv * w.cross_product(v)
	end

	def Jacobian.dq_dx_torsion(cartesians, coord_i)
		#: Arguments:
		#* cartesians - Array (4) of Coordinate
		#* coord_i - index of the coordinate (0,1,2,3)

		u = cartesians[0] - cartesians[1]
		v = cartesians[3] - cartesians[2]
		w = cartesians[2] - cartesians[1]

		lu = u.abs
		lv = v.abs
		lw = w.abs

		u.normalize!
		v.normalize!
		w.normalize!

		# Angles
		cos_u = u.dot(w)
		sin2_u = (1.0-cos_u**2)

		cos_v = v.dot(w)
		sin2_v = (1.0-cos_v**2)

		# Signs
		sign1 = Jacobian.signf(coord_i, 0, 1)
		sign2 = Jacobian.signf(coord_i, 2, 3)
		sign3 = Jacobian.signf(coord_i, 1, 2)
		d = 	sign1 / lu / sin2_u * (u.cross_product(w)) +
			sign2 / lv / sin2_v * (v.cross_product(w)) +
			sign3 * ( (u.cross_product(w)) * cos_u / lw / sin2_u -
				  (v.cross_product(w)) * cos_v / lw / sin2_v )
		return d
	end

	#=======================================================================
	# Elementary coordinates - first derivatives - numerical
	#=======================================================================

	# Numerical derivatives are included for testing purposes only

	def Jacobian.dq_dx_bond_num(cartesians, coord_i)
		#: Numerical equivalent of -->Jacobian.dq_dx_bond
		d = Coordinate.new
		3.times {|i|
			tmp = Coordinate.new
			tmp[i] = @@num_dx
			dq_plus  = cartesians[coord_i-1].distance(cartesians[coord_i]+tmp)
			dq_minus = cartesians[coord_i-1].distance(cartesians[coord_i]-tmp)
			d[i] = (dq_plus - dq_minus) / (@@num_dx * 2.0)
		}
		return d
	end

	def Jacobian.dq_dx_distinv_num(cartesians, coord_i)
		#: Numerical equivalent of -->Jacobian.dq_dx_bond
		d = Coordinate.new
		3.times {|i|
			tmp = Coordinate.new
			tmp[i] = @@num_dx
			dq_plus  = 1.0 / (cartesians[coord_i-1].distance(cartesians[coord_i]+tmp))
			dq_minus = 1.0 / (cartesians[coord_i-1].distance(cartesians[coord_i]-tmp))
			d[i] = (dq_plus - dq_minus) / (@@num_dx * 2.0)
		}
		return d
	end

	def Jacobian.dq_dx_angle_num(cartesians, coord_i)
		#: Numerical equivalent of -->Jacobian.dq_dx_angle
		d = Coordinate.new
		add = [Coordinate.new, Coordinate.new, Coordinate.new]
		3.times {|i|
			add[coord_i] = Coordinate.new
			add[coord_i][i] = @@num_dx

			dq_plus = Coordinate.angle(cartesians[0]+add[0], cartesians[1]+add[1], cartesians[2]+add[2])
			dq_minus = Coordinate.angle(cartesians[0]-add[0], cartesians[1]-add[1], cartesians[2]-add[2])
			d[i] = (dq_plus - dq_minus) / (@@num_dx * 2.0)
		}
		return d
	end

	def Jacobian.dq_dx_torsion_num(cartesians, coord_i)
		#: Numerical equivalent of -->Jacobian.dq_dx_torsion
		d = Coordinate.new
		add = [Coordinate.new, Coordinate.new, Coordinate.new, Coordinate.new]
		3.times {|i|
			add[coord_i] = Coordinate.new
			add[coord_i][i] = @@num_dx

			dq_plus = Coordinate.torsion(cartesians[0]+add[0], cartesians[1]+add[1], cartesians[2]+add[2], cartesians[3]+add[3])
			dq_minus = Coordinate.torsion(cartesians[0]-add[0], cartesians[1]-add[1], cartesians[2]-add[2], cartesians[3]-add[3])
			dq = dq_plus - dq_minus
			dq -= 2.0 * Math::PI if dq > Math::PI
			dq += 2.0 * Math::PI if dq <= -Math::PI
			d[i] = (dq) / (@@num_dx * 2.0)
		}
		return d
	end

	#=======================================================================
	# Elementary coordinates - second derivatives
	#=======================================================================
	
	def Jacobian.dq2_dxdy_bond(cartesians, coord_i, coord_j)
		#: Arguments:
		#* cartesians - Array (2) of Coordinate
		#* coord_i - index of the coordinate (0,1)
		#* coord_j - index of the coordinate (0,1)

		u = cartesians[0] - cartesians[1]
		lu = u.abs
		u.normalize!

		d = Matrix.zero(3)
		3.times{|i|
			3.times{|j|
				d[i,j] = (-1.0)**MathExt.krond_d(coord_i,coord_j) * (u[i] * u[j] - MathExt.krond_d(i,j)) / lu
			}
		}

		return d
	end

	# The formula in the paper is wrong
	#--------------------------------------------------
	# def Jacobian.dq2_dxdy_angle(cartesians, coord_i, coord_j)
	# 	#: Arguments:
	# 	#* cartesians - Array (3) of Coordinate
	# 	#* coord_i - index of the coordinate (0,1,2)
	# 	#* coord_j - index of the coordinate (0,1,2)
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	u = cartesians[0] - cartesians[1]
	# 	v = cartesians[2] - cartesians[1]
	# 	lu = u.abs
	# 	lv = v.abs
	# 	u.normalize!
	# 	v.normalize!
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	cos_a = u.dot(v)
	# 	sin_a = (1.0-cos_a**2)**0.5
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	d = Matrix.zero(3)
	# 	return d if sin_a == 0.0
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	z1 = Jacobian.signf(coord_i, 0, 1) * Jacobian.signf(coord_j, 0, 1)
	# 	z2 = Jacobian.signf(coord_i, 2, 1) * Jacobian.signf(coord_j, 2, 1)
	# 	z3 = Jacobian.signf(coord_i, 0, 1) * Jacobian.signf(coord_j, 2, 1)
	# 	z4 = Jacobian.signf(coord_i, 2, 1) * Jacobian.signf(coord_j, 0, 1)
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	3.times{|i|
	# 		3.times{|j|
	# 			d[i,j] = z1 * (u[i]*v[j] + u[j]*v[i] - 3.0*u[i]*u[j]*cos_a + MathExt.krond_d(i,j)*cos_a) / (lu**2 * sin_a) +
	# 				 z2 * (v[i]*u[j] + v[j]*u[i] - 3.0*v[i]*v[j]*cos_a + MathExt.krond_d(i,j)*cos_a) / (lv**2 * sin_a) +
	# 				 z3 * (u[i]*u[j] + v[j]*v[i] -     u[i]*v[j]*cos_a - MathExt.krond_d(i,j))       / (lu*lv * sin_a) +
	# 				 z4 * (v[i]+v[j] + u[j]*u[i] -     v[i]*u[j]*cos_a - MathExt.krond_d(i,j))       / (lu*lv * sin_a) -
	# 				 cos_a / sin_a * Jacobian.dq_dx_angle(cartesians, coord_i)[i] * Jacobian.dq_dx_angle(cartesians, coord_j)[j]
	# 		}
	# 	}
	# 	return d
	# end
	#-------------------------------------------------- 

	# Can not be implemented directly from the information in the paper
	#--------------------------------------------------
	# def Jacobian.dq2_dxdy_torsion(cartesians, coord_i, coord_j)
	# 	#: Arguments:
	# 	#* cartesians - Array (4) of Coordinate
	# 	#* coord_i - index of the coordinate (0,1,2,3)
	# 	#* coord_j - index of the coordinate (0,1,2,3)
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	d = Matrix.zero(3)
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# No coupling between terminal atoms:
	# 	if (coord_i == 0 && coord_j == 3) || (coord_i == 3 && coord_j == 0)
	# 		return d
	# 	end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	u = cartesians[0] - cartesians[1]
	# 	v = cartesians[3] - cartesians[2]
	# 	w = cartesians[2] - cartesians[1]
	# 	lu = u.abs
	# 	lv = v.abs
	# 	lw = w.abs
	# 	u.normalize!
	# 	v.normalize!
	# 	w.normalize!
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Angles
	# 	cos_u = u.dot(w)
	# 	sin_u = (1.0-cos_u**2)**0.5
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	cos_v = -v.dot(w)
	# 	sin_v = (1.0-cos_v**2)**0.5
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	3.times{|i|
	# 		3.times{|j|
	# 			d[i,j] = 0.0
	# 		}
	# 	}
	# 	return d
	# end
	#-------------------------------------------------- 

	#=======================================================================
	# Elementary coordinates - second derivatives - numerical
	#=======================================================================

	# Numerical derivatives are included for testing purposes only

	def Jacobian.dq2_dxdy_bond_num(cartesians, coord_i, coord_j)
		d = Matrix.zero(3)
		3.times{|j|
			tmp = Coordinate.new
			tmp[j] = @@num_dx 
			displaced = []
			displaced_m = []
			2.times{|x|
				if x == coord_j
					displaced << cartesians[x] + tmp
					displaced_m << cartesians[x] - tmp
				else
					displaced << cartesians[x]
					displaced_m << cartesians[x]
				end
			}
			dq_plus = Jacobian.dq_dx_bond(displaced, coord_i)
			dq_minus = Jacobian.dq_dx_bond(displaced_m, coord_i)
			v = (dq_plus - dq_minus) / (@@num_dx * 2.0)
			3.times{|i| d[i,j] = v[i]}
		}
		return d
	end

	def Jacobian.dq2_dxdy_angle_num(cartesians, coord_i, coord_j)
		d = Matrix.zero(3)
		3.times{|j|
			tmp = Coordinate.new
			tmp[j] = @@num_dx 
			displaced = []
			displaced_m = []
			3.times{|x|
				if x == coord_j
					displaced << cartesians[x] + tmp
					displaced_m << cartesians[x] - tmp
				else
					displaced << cartesians[x]
					displaced_m << cartesians[x]
				end
			}
			dq_plus = Jacobian.dq_dx_angle(displaced, coord_i)
			dq_minus = Jacobian.dq_dx_angle(displaced_m, coord_i)
			v = (dq_plus - dq_minus) / (@@num_dx * 2.0)
			3.times{|i| d[i,j] = v[i]}
		}
		return d
	end

	def Jacobian.dq2_dxdy_torsion_num(cartesians, coord_i, coord_j)
		d = Matrix.zero(3)
		3.times{|j|
			tmp = Coordinate.new
			tmp[j] = @@num_dx 
			displaced = []
			displaced_m = []
			4.times{|x|
				if x == coord_j
					displaced << cartesians[x] + tmp
					displaced_m << cartesians[x] - tmp
				else
					displaced << cartesians[x]
					displaced_m << cartesians[x]
				end
			}
			dq_plus = Jacobian.dq_dx_torsion(displaced, coord_i)
			dq_minus = Jacobian.dq_dx_torsion(displaced_m, coord_i)
			v = (dq_plus - dq_minus) / (@@num_dx * 2.0)
			3.times{|i| d[i,j] = v[i]}
		}
		return d
	end


end
