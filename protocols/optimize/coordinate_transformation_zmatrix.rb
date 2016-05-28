require "classes/internal_coordinates/internal_coordinate_set.rb"
require "protocols/optimize/coordinate_transformation.rb"
require "classes/objects/hessian_internal_coordinates.rb"

#!# Problem: optimization with dummy atoms does not work

class CoordinateTransformationZMatrix < CoordinateTransformation
	
	def initialize(geometry, settings)
		Cuby.log.puts_debug "Z-Matrix coordinate transformation"
		@z_matrix = geometry.info[:z_matrix]
		
		# Add names to all coordinates
		@z_matrix.each_coord_with_index{|coord, i|
			unless coord.name
				coord.name = "var#{i}"
			end
		}
		Cuby.log.puts_debug "Z-matrix:\n" + @z_matrix.to_s + "\n"

		@int_coords = InternalCoordinateSet.from_z_matrix(@z_matrix, geometry)

		# search for special names
		stars = 0
		excla = 0
		@int_coords.each_with_index{|ic, i|
			if ic.name =~ /\*/
				stars += 1
			elsif ic.name =~ /!/
				excla += 1
			end
		}

		if stars > 0 && excla > 0
			Cuby::error("Z-matrix can not mix coordinates labelled '*' and '!'")
		end
		
		# If there are coordinates with * in their names, save them in  @selected_coords
		# If it is nil, everything is optimized

		@selected_coords = {}
		if excla > 0 || stars == 0
			# Some coordinates labelled as not optimized - select all but those
			@int_coords.each_with_index{|ic, i|
				if ic.name !~ /!/
					if @selected_coords.has_key?(ic.name)
						# Already in the list, add this occurence
						@selected_coords[ic.name] << i
					else
						# New optimized variable, initalize array
						@selected_coords[ic.name] = [i]
					end
				end
			}
		else
			# Coordinates explicitly labelled as optimized
			@int_coords.each_with_index{|ic, i|
				if ic.name =~ /\*/
					if @selected_coords.has_key?(ic.name)
						# Already in the list, add this occurence
						@selected_coords[ic.name] << i
					else
						# New optimized variable, initalize array
						@selected_coords[ic.name] = [i]
					end
				end
			}
		end

		@dimension = @selected_coords.size
		@geometry = geometry
		@settings = settings
		
		Cuby.log.puts_debug "Optimizing #{@dimension} coordinates"
	end
	
	def geometry_to_vec(geometry)
		# Update matrices dependent on the current geometry
		@int_coords.geometry_changed
		return vector_of_selected_coordinates(@int_coords.to_vector)
	end
	
	def gradient_to_vec(gradient)
		# Convert gradient to internal coordinates
		g_int = @int_coords.b_ps_inv.transpose * gradient.to_vector
		g_int = vector_of_selected_coordinates(g_int)
		
		return g_int
	end
	
	def vec_to_geometry!(geometry, vector)
		# Update the optimized variables
		i = 0
		@selected_coords.each_pair{|name, indexes|
			indexes.each{|j|
				@z_matrix.coord_by_index(j).value = vector[i]
			}
			i += 1
		}
		
		@z_matrix.update_geometry!(geometry)
		return vector
	end
	
	def geometry_to_geometry(geometry)
		return geometry
	end
	
	def transformed_gradient(gradient)
		gi = gradient_to_vec(gradient) # Gradient in reduced internal coordinate set
		# Expand it to full dimension, filling rest with zeroes
		gfull = Vector.zero(@int_coords.size)
		i = 0
		@selected_coords.each_pair{|name, indexes|
			indexes.each{|j|
				gfull[j] = gi[i]
			}
			i += 1
		}
		
		# Convert it back to cartesian
		gcart =  @int_coords.b.transpose * @int_coords.p_matrix * gfull

		# Remove translation/rotation
		#projector = Hessian.transrot_removal_projector(@geometry)
		#gcart = projector * gcart

		return gcart
	end
	
	def initial_hessian(stat_counters = nil)   #!# not tested !!!
		unless @settings.set?(:init_hessian_lookup_file)
			return @settings[:opt_diagonal_h0]
		else
			database = HessianInternalCoordinates.from_yaml(@settings[:init_hessian_lookup_file])
			diag_hess = []
			conn = Connectivity.new(@geometry)
			all_atoms_types = conn.atom_types_simple(@geometry)
			@dimension.times do |i|
				dement = HessianInternalCoordinates::DiagonalElement.from_molecule(@geometry, self, all_atoms_types, nil, i)
				diag_hess << database.best_match_value(dement, stat_counters)
			end
			
			return Vector.from_array(diag_hess)
		end
	end

	def periodicity
		vector1 = Vector.of_size(@int_coords.size)
		vector2 = Vector.of_size(@int_coords.size)
		@int_coords.each_index{|i|
			if @int_coords[i].type == :torsion
				vector1[i] = -Math::PI
				vector2[i] = Math::PI
			else
				vector1[i] = 0.0
				vector2[i] = 0.0
			end
		}
		vector1 = vector_of_selected_coordinates(vector1)
		vector2 = vector_of_selected_coordinates(vector2)
		return [vector1, vector2]
	end
	
	#===============================================================================
	# Private
	#===============================================================================
	
	def vector_of_selected_coordinates(vec)
		subset = Vector.of_size(@dimension)
		i = 0
		@selected_coords.each_pair{|name, indexes|
			# Value is averaged (to get average gradient)
			sum = 0.0
			indexes.each{|j|
				sum += vec[j]
			}
			subset[i] = sum/indexes.size

			i += 1
		}
		return subset
	end
end
