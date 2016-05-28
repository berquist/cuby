################################################################################
#
# Module CoordinateVector
#
# Author: Jan Rezac
# Date created: 2008-08-30
# License: Cuby license
# Description: Extension of Coordinate class adding some functionality of Vector
# Status: Works, Documented
#
################################################################################

require "classes/algebra/algebra.rb"

module CoordinateVector

	#=======================================================================
	# type conversion
	#=======================================================================
	
	def to_vector # => Vector
		#: Conversion to 3D Vector
		return Vector.from_array([self.x,self.y,self.z])
	end

	#=======================================================================
	# misc
	#=======================================================================
	
	def rotate_by_matrix(rotator) # => Coordinate
		#: Rotation of vector by application of 3x3 matrix
		# Matrix - vector multiplication:
		nx = x * rotator[0,0] + y * rotator[0,1] + z * rotator[0,2]
		ny = x * rotator[1,0] + y * rotator[1,1] + z * rotator[1,2]
		nz = x * rotator[2,0] + y * rotator[2,1] + z * rotator[2,2]
		return Coordinate[nx,ny,nz]
	end

	def rotate_by_matrix!(rotator) # => nil
		#: In-place rotation of vector by application of 3x3 matrix
		set_coord(self.rotate_by_matrix(rotator))
		return nil
	end

	def CoordinateVector.rotator_from_angle_axis(angle,axis) # => Matrix
		#: Creates rotator matrix from ange (float) and axis (3D vector)
		return Algebra::Vector3DMethods.rotator_from_angle_axis(angle,axis.to_vector)
	end

	def rotate_angle_axis(angle,axis) # => Coordinate
		#: Returns vector rotated by angle (float) around axis (3D vector)
		r = CoordinateVector.rotator_from_angle_axis(angle,axis)
		return rotate_by_matrix(r)
	end

	def rotate_angle_axis!(angle,axis) # => nil
		#: Rotates vector by angle (float) around axis (3D vector)
		r = CoordinateVector.rotator_from_angle_axis(angle,axis)
		rotate_by_matrix!(r)
		return nil
	end
end

#===============================================================================
# Extend 3D vector to allow conversion to coordinate
#===============================================================================

module Vector3DMethods
	def to_coordinate # => Coordinate
		#: Conversion of 3D vector to coordinate
		raise_if_not_3D
		return Coordinate.new(x,y,z)
	end
end
