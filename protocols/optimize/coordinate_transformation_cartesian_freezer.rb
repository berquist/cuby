require "protocols/optimize/coordinate_transformation.rb"

class CoordinateTransformationCartesianFreezer < CoordinateTransformation
	def initialize(geometry, settings)
		@settings = settings
		Cuby.log.puts_debug "Caretesian coordinate freezer active"
		# Build list of optimized coordinates
		frozen_coords = []
		@settings[:opt_freeze_cartesian].each_pair{|index, coords|
			coords = coords.strip.downcase
			frozen_coords << (index-1)*3 + 0 if coords =~ /x/
			frozen_coords << (index-1)*3 + 1 if coords =~ /y/
			frozen_coords << (index-1)*3 + 2 if coords =~ /z/
		}
		@optimized_coords = (0..(geometry.size*3-1)).to_a - frozen_coords
		Cuby.log.puts_debug "#{frozen_coords.size} coordinates frozen"

		# Save values of frozen coordinates
		# Not needed in the current implementation
		#--------------------------------------------------
		# full_vec = geometry.to_vector
		# @frozen_values = {}
		# frozen_coords.each{|i|
		# 	@frozen_values[i] = full_vec[i]
		# }
		#-------------------------------------------------- 

		@dimension = @optimized_coords.size
	end

	def geometry_to_vec(geometry)
		# Vector of values of optimized coordinates
		full_vec = geometry.to_vector
		a = []
		@optimized_coords.each{|i|
			a << full_vec[i]
		}
		return Vector.from_array(a)
	end

	def gradient_to_vec(gradient)
		# Vector of gradients of optimized coordinates
		full_vec = gradient.to_vector
		a = []
		@optimized_coords.each{|i|
			a << full_vec[i]
		}
		return Vector.from_array(a)
	end

	def vec_to_geometry!(geometry, vector)
		# Update the optimized coordinates in geometry
		@optimized_coords.each_index{|i|
			a_i = @optimized_coords[i] / 3
			c_i = @optimized_coords[i] % 3
			geometry[a_i][c_i] = vector[i]
		}
		return vector
	end

	def geometry_to_geometry(geometry)
		return geometry
	end
	
	def transformed_gradient(gradient)
		# Zero-out the gradient on non-moving atoms
		gradient2 = Vector.zero(gradient.size)
		@optimized_coords.each{|i|
			gradient2[i] = gradient[i]
		}
		return gradient2
	end
	
	def hessian_to_matrix(hessian)
		# Zero-out the Hessian on non-moving atoms
		hessian2 = Matrix.zero(hessian.m, hessian.n)
		@optimized_coords.each{|i|
			@optimized_coords.each{|j|
				hessian2[i,j] = hessian[i,j]
			}
		}
		return hessian2
	end
	
	def initial_hessian(stat_counters = nil)
		return @settings[:opt_diagonal_h0]
	end
end
