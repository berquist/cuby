require "protocols/optimize/coordinate_transformation.rb"

class CoordinateTransformationNone < CoordinateTransformation
	def initialize(geometry, settings)
		Cuby.log.puts_debug "No coordinate transformation - optimization in complete cartesian space"
		@dimension = geometry.size * 3
		@settings = settings
	end

	def geometry_to_vec(geometry)
		return geometry.to_vector
	end

	def gradient_to_vec(gradient)
		return gradient.to_vector
	end

	def vec_to_geometry!(geometry, vector)
		geometry.update_from_vector!(vector)
		return vector
	end

	def geometry_to_geometry(geometry)
		return geometry
	end
	
	def transformed_gradient(gradient)
		return gradient
	end
	
	def hessian_to_matrix(hessian)
		return hessian
	end
	
	def initial_hessian(stat_counters = nil)
		return @settings[:opt_diagonal_h0]
	end
end
