require "classes/internal_coordinates/internal_coordinate_set.rb"
require "protocols/optimize/coordinate_transformation_redundant.rb"
require "classes/objects/hessian_internal_coordinates.rb"

class CoordinateTransformationRedundantFreezer < CoordinateTransformationRedundant

	def initialize(geometry, settings)
		@atomlist = geometry.atomlist_from_selection(settings[:optimize_region])
		super(geometry, settings)
		
		Cuby.log.puts_debug "Freezer in redundant coordinates."
	end
	
	def gradient_to_vec(gradient)
		# Convert cartesian gradient to the optimization coordinates
		# Cartesian to internal:
		g_trunc = []
		g_vec = gradient.to_vector
		@atomlist.each{|ai|
			3.times{|k|
				g_trunc << g_vec[3*ai+k]
			}
		}
		if @settings[:development][:sparse_mode]
			g_int = @int_coords.csr_b * @int_coords.csr_btb.solve(Vector.from_array(g_trunc))
		else
			g_int = @int_coords.b_ps_inv.transpose * Vector.from_array(g_trunc)
		end
		

		# Done when full space is optimized
		return g_int if @int_coords.size == @coords_list.size
		
		# Extract optimized subspace
		vector = Vector.of_size(@dimension)
		@coords_list.each_index{|i|
			vector[i] = g_int[@coords_list[i]]
		}
		return vector
	end

	def vec_to_geometry!(geometry, vector)
		# Reconstruct coordinates in full space 
		if @int_coords.size != @coords_list.size
			# fill vector
			a = []
			j = 0
			@int_coords.each_index{|i|
				if @coords_list.include?(i)
					a[i] = vector[j]
					j += 1
				else
					a[i] = @starting_values[i]
				end
			}
			vector = Vector.from_array(a)
		end
		
		# Iterative transformation from internal to cartesian
		q0 = @int_coords.to_vector
		sq = vector - q0
		fix_torsions!(sq, @int_coords)
		
		sqdbg = sq * 1.0   #!# dbg
		
		if @settings[:development][:sparse_mode]   # Projection onto physically valid step
			sq = @int_coords.csr_b * @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
		else
			sq = @int_coords.p_matrix * sq
		end
		
		#~ unless sq.abs == 0.0 or sqdbg.abs == 0.0
			#~ Cuby::log.puts_debug("Anglength sq0--projection: #{Math::acos(sqdbg.dot(sq)/sq.abs/sqdbg.abs)*180.0/Math::PI}deg, (#{sqdbg.abs}, #{sq.abs})")
		#~ end
		
		sq *= 0.2/sq.max_abs unless sq.max_abs == 0.0
		
		sq0_init = sq * 1.0
		sq0 = sq * 1.0
		vector = q0 +  sq0_init
		
		x0 = geometry.to_vector
		x0_trunc = []   # slice cartesian vector by atomlist
		@atomlist.each{|ai|
			3.times{|k|
				x0_trunc << x0[3*ai+k]
			}
		}
		x0 = Vector.from_array(x0_trunc)
		x00 = Vector.from_array(x0_trunc)
		
		if @settings[:development][:sparse_mode]
			x1_init = x0 + @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
		else
			x1_init = x0 + @int_coords.b_ps_inv * sq
		end
		
		upd_x1_init = geometry.to_vector
		@atomlist.each_index{|i|
			ai = @atomlist[i]
			3.times{|k|
				upd_x1_init[3*ai+k] = x1_init[3*i+k]
			}
		}
		i = 0
		while i < 25 do
			if @settings[:development][:sparse_mode]
				x1 = x0 + @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
			else
				x1 = x0 + @int_coords.b_ps_inv * sq
			end
			
			upd_x = geometry.to_vector
			@atomlist.each_index{|i|
				ai = @atomlist[i]
				3.times{|k|
					upd_x[3*ai+k] = x1[3*i+k]
				}
			}
			geometry.update_from_vector!(upd_x)
			if (x1 - x0).abs < 1.0e-8   #!# tune? dependency on stepsize?
				@int_coords.geometry_changed
				Cuby::log.puts_debug "Transformation: dx.abs < 1.0e-8, iterations=#{i}; |sq| = #{sq.abs}, |dx| = #{(x1-x0).abs}"
				break
			end
			q1 = @int_coords.to_vector
			sq = vector - q1
			fix_torsions!(sq, @int_coords)
			if sq.abs > sq0_init.abs
				@int_coords.geometry_changed
				Cuby::log.puts_debug "Geometry changed: sq > sq_init !!!"
				if @settings[:development][:sparse_mode]   # Projection onto physically valid step
					sq = @int_coords.csr_b * @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
				else
					sq = @int_coords.p_matrix * sq
				end
				i += 3   # iteration is poor, cannot do the same all 25 steps
			end
			if sq.abs < 1.0e-12   #!# tune? dependency on stepsize?
				Cuby::log.puts_debug "Transformation: dq.abs < 1.0e-12, iteration=#{i}; |sq| = #{sq.abs}, |dx| = #{(x1-x0).abs}"
				break
			end
			i += 1
			
			#~ unless sq.abs == 0.0 or sq0.abs == 0.0
				#~ Cuby::log.puts_debug("Anglength sq0--sq (#{i}. step): #{Math::acos(sq0.dot(sq)/sq.abs/sq0.abs)*180.0/Math::PI}deg, (#{sq0.abs}, #{sq.abs})")
			#~ end
			Cuby::log.puts_debug "Transformation iterations reached #{i} !!! |sq| = #{sq.abs}, |dx| = #{(x1-x0).abs}" if i > 24
			
			x0 = x1
			sq0 = sq
			if i % 5 == 0   # update B matrix every 5 steps at least
				@int_coords.geometry_changed
				Cuby::log.puts_debug "Geometry changed: 5th step."
				if @settings[:development][:sparse_mode]   # Projection onto physically valid step
					sq = @int_coords.csr_b * @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
				else
					sq = @int_coords.p_matrix * sq
				end
			end
		end
		
		Cuby::log.puts_debug "Step cartesian: #{(x1-x00).max_abs}"
		return @int_coords.to_vector
	end
	
	def transformed_gradient(gradient)
		# Build cartesian gradient corresponding to the one used in
		# optimization. This gradient is used to check the convergence
		g_trunc = []
		g_vec = gradient.to_vector
		@atomlist.each{|ai|
			3.times{|k|
				g_trunc << g_vec[3*ai+k]
			}
		}
		g_trunc_vec = Vector.from_array(g_trunc)
		
		# No freezing -> return gradient as is
		return g_trunc_vec if @int_coords.size == @coords_list.size
		
		# Otherwise, transform it to internal coords
		if @settings[:development][:sparse_mode]
			g_int = @int_coords.csr_b * @int_coords.csr_btb.solve(g_trunc_vec)
		else
			g_int = @int_coords.b_ps_inv.transpose * g_trunc_vec
		end
		
		g_int.each_index{|i|
			g_int[i] = 0.0 unless @coords_list.include?(i)
		}
		# Return gradient transformed back to cartesian coordinates
		if @settings[:development][:sparse_mode]
			# as soon as B matrix has independent columns, B^T P g = B^T B B^+ g = (B^T B) (B^T B)^(-1) B^T g = B^T g
			return @int_coords.csr_bt * g_int
		else
			return @int_coords.b.transpose * g_int
		end
		
	end
	
end
