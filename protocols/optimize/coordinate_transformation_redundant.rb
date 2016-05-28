require "classes/internal_coordinates/internal_coordinate_set.rb"
require "protocols/optimize/coordinate_transformation.rb"
require "protocols/optimize/hessian_estimate_redundant.rb"
require "classes/objects/hessian_internal_coordinates.rb"
require "classes/objects/hessian.rb"

class CoordinateTransformationRedundant < CoordinateTransformation

	attr_reader :int_coords
	
	def initialize(geometry, settings)
		# Initialization of the coordinate transformation
		
		Cuby.log.puts_debug "Redundant coordinate transformation"
		# Build a complete set of redundant coordinates
		@int_coords = InternalCoordinateSet.redundant(geometry, settings)
		# Make a lsit of coordinates that are optimized
		@coords_list = build_coords_list(settings)
		@dimension = @coords_list.size
		@settings = settings
		@geometry = geometry
		# Save initial values of the coordinates
		@starting_values = []
		@int_coords.each_index{|i|
			@starting_values[i] = @int_coords[i].value
		}
		# Print information on internal coordinates
		Cuby::log.puts_v(:normal, "Internal coordinates used:")
		@int_coords.count_type_name.each_pair{|type, count|
			Cuby::log.puts_v(:normal, sprintf("%8d   %s", count, type.to_s))
		}
		Cuby::log.puts_v(:normal, "")

		if @settings[:development][:opt_plot_matrices]
			puts "B matrix"
			@int_coords.b.plot_gnuplot('matrix', true, "eog")
			puts "B matrix pseudoinverse"
			@int_coords.b_ps_inv.plot_gnuplot('matrix', true, "eog")
			puts "B.trans * B"
			m2 = @int_coords.b.transpose * @int_coords.b
			m2.plot_gnuplot('matrix', true, "eog")
			puts "Sizes:"
			puts "B: #{@int_coords.b.m * @int_coords.b.n}"
			puts "B.trans * B: #{m2.m * m2.n}"
			nz = m2.count_nonzero
			puts "out of it, nonzero: #{nz} (ratio: #{nz.to_f / (m2.m * m2.n)})"
			exit
		end
	end
	
	def geometry_to_vec(geometry)
		# Convert cartesian coordinates to the optimization coordinates
		# stored as a vector

		# Cartesian to internal, done when full space is optimized
		return @int_coords.to_vector if @int_coords.size == @coords_list.size
		
		# Extract optimized subspace
		vector = Vector.of_size(@dimension)
		@coords_list.each_index{|i|
			vector[i] = @int_coords[@coords_list[i]].value
		}
		return vector
	end

	def gradient_to_vec(gradient)
		# Convert cartesian gradient to the optimization coordinates
		# Cartesian to internal:
		#!# TODO
		#~ if @settings[:development][:sparse_mode]
			#~ 
		#~ else
			g_int = @int_coords.b_ps_inv.transpose * gradient.to_vector
		#~ end
		
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
		# Translation/rotation removal - not used
		# @transrot_projector = Hessian.transrot_removal_projector(geometry)
		@transrot_projector = 1.0

		# Reconstruct coordinates in full space 
		if @int_coords.size != @coords_list.size
			# doplnit vector
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

		#!# TODO
		#~ if @settings[:development][:sparse_mode]   # Projection onto physically valid step
			#~ 
		#~ else
			sq = @int_coords.p_matrix * sq
		#~ end

		sq0 = sq0_init = sq
		vector = q0 +  sq0_init
		
		x0 = geometry.to_vector
		#!# TODO
		#~ if @settings[:development][:sparse_mode]
			#~ x1_init = x0 + @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
		#~ else
			x1_init = x0 + @transrot_projector * @int_coords.b_ps_inv * sq
		#~ end
		
		upd_x1_init = x1_init
		
		i = 0
		while i < 25 do
			#!# TODO
			#~ if @settings[:development][:sparse_mode]
				#~ x1 = x0 + @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
			#~ else
				x1 = x0 + @transrot_projector * (@int_coords.b_ps_inv * sq)
			#~ end
			upd_x = x1
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
				#!# TODO
				#~ if @settings[:development][:sparse_mode]   # Projection onto physically valid step
					#~ sq = @int_coords.csr_b * @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
				#~ else
					sq = @int_coords.p_matrix * sq
				#~ end
			end
			if sq.abs < 1.0e-12   #!# tune? dependency on stepsize?
				Cuby::log.puts_debug "Transformation: dq.abs < 1.0e-12, iteration=#{i}; |sq| = #{sq.abs}, |dx| = #{(x1-x0).abs}"
				break
			end
			i += 1
			Cuby::log.puts_debug "Transformation iterations reached #{i} !!! |sq| = #{sq.abs}, |dx| = #{(x1-x0).abs}" if i > 24
			x0 = x1
			sq0 = sq
			if i % 5 == 0   # update B matrix every 5 steps at least
				@int_coords.geometry_changed
				#!# TODO
				#~ if @settings[:development][:sparse_mode]   # Projection onto physically valid step
					#~ sq = @int_coords.csr_b * @int_coords.csr_btb.solve(@int_coords.csr_bt * sq)
				#~ else
					sq = @int_coords.p_matrix * sq
				#~ end
			end
		end
		
		return @int_coords.to_vector
	end

	def geometry_to_geometry(geometry)
		# Complete geometry is optimized, no changes has to be made
		return geometry
	end
	
	def transformed_gradient(gradient)
		# Build cartesian gradient corresponding to the one used in
		# optimization. This gradient is used to check the convergence

		# No freezing -> return gradient as is
		return gradient if @int_coords.size == @coords_list.size

		# Otherwise, transform it to internal coords
		#!# TODO
		#~ if @settings[:development][:sparse_mode]
			#~ g_int = @int_coords.csr_b * @int_coords.csr_btb.solve(g_trunc_vec)
		#~ else
			g_int = @int_coords.b_ps_inv.transpose * gradient.to_vector
		#~ end
		
		g_int.each_index{|i|
			g_int[i] = 0.0 unless @coords_list.include?(i)
		}
		
		# Return gradient transformed back to cartesian coordinates
		#!# TODO
		#~ if @settings[:development][:sparse_mode]
			#~ # as soon as B matrix has independent columns, B^T P g = B^T B B^+ g = (B^T B) (B^T B)^(-1) B^T g = B^T g
			#~ return @int_coords.csr_bt * g_int
		#~ else
			return @transrot_projector * @int_coords.b.transpose * @int_coords.p_matrix * g_int
		#~ end
	end
	
	def hessian_to_matrix(hessian, gradient)
		return @int_coords.b_ps_inv.transpose * (hessian - @int_coords.k_matrix(gradient)) * @int_coords.b_ps_inv
		# symbolically: H_q = B_ps_inv^T (H_x - K) B_ps_inv
	end
	
	def initial_hessian(stat_counters = nil)
		case @settings[:hessian_estimate]
		when :single_number
			# The simplest estimate: unit matrix * some number
			return @settings[:opt_diagonal_h0]
		when :by_coordinate
			return HessianEstimateRedundant.simple(@int_coords, @coords_list)
		when :from_lookup
			# Vector of elements from an atom-type based precalculated lookup
			database = HessianInternalCoordinates.from_yaml(@settings[:init_hessian_lookup_file])
			diag_hess = []
			conn = Connectivity.new(@geometry)
			all_atoms_types = conn.atom_types_simple(@geometry)
			@coords_list.each do |i|
				dement = HessianInternalCoordinates::DiagonalElement.from_molecule(@geometry, self, all_atoms_types, nil, i)
				diag_hess << database.best_match_value(dement, stat_counters)
			end
			return Vector.from_array(diag_hess)
		when :fischer_almlof
			return HessianEstimateRedundant.fischer_almlof(@int_coords, @coords_list, @geometry)
		when :lindh
			return HessianEstimateRedundant.lindh(@int_coords, @coords_list, @geometry)
		when :swart
			return HessianEstimateRedundant.swart(@int_coords, @coords_list, @geometry)
		when :from_file
			# Does not work with freezing atoms
			if @settings[:optimize_region] and @settings[:optimize_region] != "%all()"
				Cuby::error("Cartesian hessian does not work with :from_file and :optimize_region.")
			end
			# Load a hessian file
			# Convert it to redundant coordinates (assume gradient = 0)
			f = File.open(@settings[:hessian_read], "r")
			cart_hess = YAML.load(f)
			f.close
			return hessian_to_matrix(cart_hess, Vector.zero(@geometry.size*3))
		when :from_file_shift
			# Does not work with freezing atoms
			if @settings[:optimize_region] and @settings[:optimize_region] != "%all()"
				Cuby::error("Cartesian hessian does not work with :from_file and :optimize_region.")
			end
			# Load a hessian file
			# Convert it to redundant coordinates (assume gradient = 0) and make it non-singular with min_hess value
			f = File.open(@settings[:hessian_read], "r")
			cart_hess = YAML.load(f)
			f.close
			
			int_hess = hessian_to_matrix(cart_hess, Vector.zero(@geometry.size*3))
			eig_sys = int_hess.eigensystem
			
			min_hess = 1.0
			min = min_hess
			eig_sys[1] = eig_sys[1].to_vector
			eig_sys[1].each{|l| min = l if l < min}
			if min < min_hess
				eig_sys[1].each_index{|i| eig_sys[1][i] += min_hess - min}
				return eig_sys[0] * Matrix.diagonal(eig_sys[1]) * eig_sys[0].inverse
			else
				return int_hess
			end
		end
	end

	def periodicity
		vector1 = Vector.of_size(@dimension)
		vector2 = Vector.of_size(@dimension)
		@coords_list.each_index{|i|
			if @int_coords[@coords_list[i]].type == :torsion
				vector1[i] = -Math::PI
				vector2[i] = Math::PI
			else
				vector1[i] = 0.0
				vector2[i] = 0.0
			end
		}
		return [vector1, vector2]
	end
	
	#=======================================================================
	# Private methods
	#=======================================================================
	
	def build_coords_list(settings)
		# Make a list of coordinates that are optimized
		# uses following keywords:
		# * redcoord_freeze_class
		cl = []
		@int_coords.each_index{|i|
			cl << i unless
				settings[:redcoord_freeze_class].include?(:bonds) && @int_coords[i].type == :distance ||
				settings[:redcoord_freeze_class].include?(:angles) && @int_coords[i].type == :angle ||
				settings[:redcoord_freeze_class].include?(:torsions) && @int_coords[i].type == :torsion
		}
		if cl.empty?
			Cuby::error("There is nothing to be optimized!\nTry to remove some coordinate-type from \"redcoord_freeze_class\"")
		end
		return cl
	end
	
	def fix_torsions!(q, ic)
		q.each_index{|i|
			if ic[i].type == :torsion
				q[i] = q[i] % (Math::PI * 2) if q[i] > Math::PI * 2
				q[i] -= 2.0*Math::PI if q[i] > Math::PI

				q[i] = -((-q[i]) % (Math::PI * 2)) if q[i] < -Math::PI * 2
				q[i] += 2.0*Math::PI if q[i] <= -Math::PI
			end
		}
	end
	
end
