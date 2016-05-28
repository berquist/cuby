require "protocols/optimize/coordinate_transformation.rb"

class CoordinateTransformationFreezer < CoordinateTransformation
	def initialize(geometry, settings)
		Cuby.log.puts_debug "Cartesian freezer - optimization in reduced cartesian space"
		@atomlist = geometry.atomlist_from_selection(settings[:optimize_region])
		@settings = settings
		@dimension = @atomlist.size * 3
	end

	def geometry_to_vec(geometry)
		a = []
		@atomlist.each{|i|
			geometry[i].each{|c| a << c}
		}
		return Vector.from_array(a)
	end

	def gradient_to_vec(gradient)
		a = []
		@atomlist.each{|i|
			gradient[i].each{|c| a << c}
		}
		return Vector.from_array(a)
	end

	def vec_to_geometry!(geometry, vector)
		@atomlist.each_index{|i|
			geometry[@atomlist[i]].x = vector[i*3]
			geometry[@atomlist[i]].y = vector[i*3+1]
			geometry[@atomlist[i]].z = vector[i*3+2]
		}
		return vector
	end

	def geometry_to_geometry(geometry)
		return geometry.geometry_from_list(@atomlist)
	end
	
	def transformed_gradient(gradient)
		newgrad = Vector.of_size(@dimension)
		@atomlist.each_index{|i|
			newgrad[i*3]	 = gradient[@atomlist[i]*3]
			newgrad[i*3 + 1] = gradient[@atomlist[i]*3 + 1]
			newgrad[i*3 + 2] = gradient[@atomlist[i]*3 + 2]
		}
		return newgrad
	end
	
	def hessian_to_matrix(hessian)
		h = Matrix.zero(3*@atomlist.size,3*@atomlist.size)
		@atomlist.each_index{|i|
		3.times{|k|
			@atomlist.each_index{|j|
			3.times{|l|
				h[3*i+k,3*j+l] = hessian(3*@atomlist[i]+k,3*@atomlist[j]+l)
			}}
		}}
	end
	
	def initial_hessian(stat_counters = nil)
		return @settings[:opt_diagonal_h0]
	end
end
