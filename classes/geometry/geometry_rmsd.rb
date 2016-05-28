################################################################################
#
# Module GeometryRMSD
#
# Author: Jan Rezac
# Date created: 2008-11-04
# License: Cuby license
# Description: Methods for RMSD calculations and fitting geometries
# Extends: Geometry
# Status: Works, Documented
#
################################################################################

#: Calculation of RMSD between two geometries and quaternion fit to minimize it.

module GeometryRmsd

	#=======================================================================
	# RMSD calculation 
	#=======================================================================
	
	def rmsd_noshift(geometry) # => Float
		#: RMSD calculation withou translation
		sum = 0.0
		size.times{|i|
			3.times {|j|
				sum += (self[i][j] - geometry[i][j])**2
			}
		}
		sum = (sum / size )**0.5
		return sum
	end

	def rmsd(geometry) # => Float
		#: Calculation of RMSD after translation to centers of mass of the geometries
		center_this  = weighted_center(:unit)
		center_other = geometry.weighted_center(:unit)

		shifted_this  = translate(center_this * -1.0)
		shifted_other = geometry.translate(center_other * -1.0)
		
		return shifted_this.rmsd_noshift(shifted_other)
	end

	def fit_rmsd(template_geo) # => Geometry
		#: Fit to minimize RMSD by translation and rotation of the geometry
		center_this  = weighted_center(:unit)
		center_template = template_geo.weighted_center(:unit)

		shifted_this  = translate(center_this * -1.0)
		shifted_template = template_geo.translate(center_template * -1.0)

		result = shifted_this.quatfit(shifted_template)

		result.translate!(center_template)

		return result
	end

	#===============================================================================
	# Quaternion fit for minimizing the RMSD
	#===============================================================================

	def quatfit_rotator(template)
		#: Quaternion fit
		#: Based on Karney C.F.F, Journal of molecular graphics and modelling 25 (2007), 595-604
		geo1 = template
		geo2 = self

		# Construct the equation to fit
		b_mat = Matrix.zero(4,4)
		geo1.each_index{|i|
			# Geo = y
			# geo2 = x
			a = geo1[i] + geo2[i]
			b = geo1[i] - geo2[i]

			a_mat = Matrix[
				[0,	-b.x,	-b.y,	-b.z],
				[b.x,	0,	-a.z,	a.y],
				[b.y,	a.z,	0,	-a.x],
				[b.z,	-a.y,	a.x,	0]
			]

			b_mat = b_mat + a_mat.transpose * a_mat
		}
		b_mat = b_mat * (1.0/geo1.size)

		# Find lowest eigenvalue
		eigenvec, eigenval, eigenval_i =  b_mat.eigensystem
		i = 0
		mineig = eigenval[0,0]
		4.times { |j|
			if eigenval[j,0] < mineig
				i = j
				mineig = eigenval[j,0]
			end
		}

		# Lowest eigenvector is the quaternion of the desired rotation
		minvec = eigenvec.column_as_array(i)
		q0, q1, q2, q3 = minvec

		# Convert quaternion representation of the rotation to rotation matrix
		r_mat = Matrix[	[1.0 - 2.0*q2**2 - 2.0*q3**2,	2.0*q1*q2 - 2.0*q0*q3,		2.0*q1*q3 + 2.0*q0*q2],
		       		[2.0*q2*q1 + 2.0*q0*q3,		1.0 - 2.0*q3**2 - 2.0*q1**2,	2.0*q2*q3 - 2.0*q0*q1],
		       		[2.0*q3*q1 - 2.0*q0*q2,		2.0*q3*q2 + 2.0*q0*q1,		1.0 - 2.0*q1**2 - 2.0*q2**2]]

		return r_mat
	end

	def quatfit(template)
		# Create new geometry as a copy of the original one,
		# rotate the atoms
		r_mat = quatfit_rotator(template)
		result = self.deep_copy
		result.each {|atom| atom.rotate_by_matrix!(r_mat)}

		return result
	end
end	
