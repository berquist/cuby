#: The CoordinateTransformation objects provide translation between the full cartesian
#: coordinates use for the calculation and the coordinate space being optimized.

class CoordinateTransformation
	attr_reader :dimension

	def initialize(geometry, settings)
		# Initialization of the coordinate transformation
		@settings = settings
		raise "This is just a template, the method has to be redefined in descendant class"
	end

	def geometry_to_vec(geometry)
		# Convert cartesian coordinates to the optimization coordinates
		# stored as a vector
		raise "This is just a template, the method has to be redefined in descendant class"
	end

	def gradient_to_vec(gradient)
		# Convert cartesian gradient to the optimization coordinates
		raise "This is just a template, the method has to be redefined in descendant class"
	end

	def vec_to_geometry!(geometry, vector)
		# Update geometry with coordinates from a vector in the 
		# optimization coordinate system
		raise "This is just a template, the method has to be redefined in descendant class"
		# The vector might be altered by the transformation, return the actual one
		return vector
	end
	
	def geometry_to_geometry(geometry)
		# Return a subset of the geometry corresponding to the
		# optimized geometry (used by transformations that freeze
		# cartesian coordinates by not including them in the geometry)
		raise "This is just a template, the method has to be redefined in descendant class"
	end

	def transformed_gradient(gradient)
		# Returns cartesian gradient corresponding to the one which
		# is used in the optimization. This gradient is used for
		# evaluation of the convergence.
		raise "This is just a template, the method has to be redefined in descendant class"
	end
	
	def hessian_to_matrix(hessian, gradient = nil)
		# Convert cartesian hessian to the optimization coordinates
		# Redundant coordinates need also gradient
		raise "This is just a template, the method has to be redefined in descendant class"
	end
	
	def initial_hessian(stat_counters = nil)
		# Returns numeric type as dimension or vector as estimate of Hessian diagonal
		raise  "This is just a template, the method has to be redefined in descendant class"
	end

	def periodicity
		# Periodicity, for periodic coordinates it should yield a vector
		return nil
	end
end
