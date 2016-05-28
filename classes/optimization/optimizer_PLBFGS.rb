require "classes/optimization/optimizer.rb"
require "classes/optimization/optimizer_common.rb"

require "set"
class SparseMatrix

	def color_array
		ca = []
		@m.times{|i|
			tmp = Set.new
			row_each_nonzero_index(i){|j|
				row_each_nonzero_index(j){|k|
					tmp.add(ca[k]) if k != i && ca[k]
				}
				tmp.add(ca[j]) if ca[j]
			}
			c = 0
			c += 1 while tmp.include?(c)
			ca[i] = c
		}
		return ca
	end

end

class OptimizerPLBFGS < Optimizer
	include OptimizerCommon
	
	
	def get_geometry(geometry)
		@n_atoms = geometry.size
		@conn = Connectivity.new(geometry)

		# H-bonds
		#--------------------------------------------------
		# @conn.h_bond_list(geometry).each{|hba|
		# 	@conn.add_bond(hba[1], hba[2])
		# }
		#-------------------------------------------------- 

		@ca = @conn.color_array

		@n_colors = 0
		@ca.each_index{|i|
			@n_colors = @ca[i] if @n_colors < @ca[i]
		}
		@n_colors += 1
		@stat_counters[:colors] = @n_colors

		if @settings[:plbfgs_colors] == -1
			@settings[:plbfgs_colors] = @n_colors
		end

		# Initialize hessian
		@h0 = SparseMatrix.new(geometry.size * 3)
	end
	
	def initialize(settings)
		# Generic initialization
		super(settings)

		# Trust radius initialization
		@maxstep = @settings[:trust_radius]
		@maxstep_max = @settings[:trust_radius_max]
		@maxstep_min = @settings[:trust_radius_min]

		# l-BFGS setup
		@lbfgs_n = @settings[:lbfgs_n]
		
		@oldgrad = []
		@oldvec = []
		
		@step_scaling = @settings[:opt_tr_mode]
		@step_scaling = :abs if @step_scaling == :default
		
		@idx = 0   # index in color list
		@dif_grad = []
		@diff = @settings[:plbfgs_init_d]   # d for finite differences
	end
	
	def parallel
		return @n_colors * 3 + 1
	end
	
	def direction_queue
		# Build and queue calculations
		if @cycle == 0
			# In first cycle, all colors are calculated
			@n_colors.times{|color|
				3.times{|xyz|
					c_vec = Vector.zero(@vector.size)
					@n_atoms.times{|j|
						c_vec[3*j+xyz] = @diff if @ca[j] == color
					}
					i = 3*color+xyz
					@calculation_queue_proc.call(i+1,@vector + c_vec)
				}
			}
		else
			# In subsequent cycles, only selected no. of colors is updated
			@settings[:plbfgs_colors].times{|ii|
				color = (@idx + ii) % @n_colors
				
				3.times{|xyz|
					c_vec = Vector.zero(@vector.size)
					@n_atoms.times{|j|
						c_vec[3*j+xyz] = @diff if @ca[j] == color
					}
					i = 3*color+xyz
					@calculation_queue_proc.call(i+1,@vector + c_vec)
				}
			}
		end
	end
	
	def direction(oldgrad, oldvec)
		# record size
		m = oldgrad.size - 1
		m = @lbfgs_n if m > @lbfgs_n
		
		# LBFGS step calculation
		q = oldgrad[0]
		s = []
		y = []
		ro = []
		alpha = []
		i = 1
		while i < m # Iterate over all saved records
			s[i] = oldvec[i-1] - oldvec[i]
			y[i] = oldgrad[i-1] - oldgrad[i]
			ro[i] = 1.0 / (y[i].dot(s[i]))
			alpha[i] = ro[i] * (s[i].dot(q))

			q -= alpha[i] * y[i]

			i += 1
		end
		
		if @cycle == 0
			# Whole hessian is updated
			@n_colors.times{|color|
				# Calculate differences
				3.times{|xyz|
					i = 3*color+xyz
					# Use previously queued calculations:
					@dif_grad[3*color+xyz] = @calculation_run_proc.call(i+1)[1] - @gradient   # moves "color" atoms in direction "xyz"
					@stat_counters[:energy_gradient_calls] += 1
				}
			}
			@n_colors.times{|color|
				# Update hessian with calculated differences
				@n_atoms.times{|j|
					if @ca[j] == color   # case 1: color of atom "j" is "color"
						3.times{|k| 3.times{|l|
							# fills 3x3 diagonal block
							@h0[3*j+k,3*j+l] = @h0[3*j+l,3*j+k] = (@dif_grad[3*color+l][3*j+k]) / @diff
						}}
					else   # case 2: color of atom "j" is different from "color"
						@conn.bound_atoms(j).each{|i|
							if @ca[i] == color   # if it finds neighbour atom with color "color", fills corresponding 3x3 block
								3.times{|k| 3.times{|l|
									@h0[3*i+k,3*j+l] = @h0[3*j+l,3*i+k] = (@dif_grad[3*color+k][3*j+l]) / @diff
								}}
								break
							end
						}
					end
				}
			}
		else
			# Selected number of colors is updated
			@settings[:plbfgs_colors].times{|ii|
				color = (@idx + ii) % @n_colors
				# Calculate differences
				3.times{|xyz|
					i = 3*color+xyz
					# Use previously queued calculations:
					@dif_grad[3*color+xyz] = @calculation_run_proc.call(i+1)[1] - @gradient   # moves "color" atoms in direction "xyz"
					@stat_counters[:energy_gradient_calls] += 1
				}
			}
			@settings[:plbfgs_colors].times{|ii|
				color = (@idx + ii) % @n_colors
				# Update hessian with calculated differences
				@n_atoms.times{|j|
					if @ca[j] == color   # case 1: color of atom "j" is "color"
						3.times{|k| 3.times{|l|
							# fills 3x3 diagonal block
							@h0[3*j+k,3*j+l] = @h0[3*j+l,3*j+k] = (@dif_grad[3*color+l][3*j+k]) / @diff
						}}
					else   # case 2: color of atom "j" is different from "color"
						@conn.bound_atoms(j).each{|i|
							if @ca[i] == color   # if it finds neighbour atom with color "color", fills corresponding 3x3 block
								3.times{|k| 3.times{|l|
									@h0[3*i+k,3*j+l] = @h0[3*j+l,3*i+k] = (@dif_grad[3*color+k][3*j+l]) / @diff
								}}
								break
							end
						}
					end
				}
			}
			
			@idx += @settings[:plbfgs_colors]
			@idx %= @n_colors
		
		end
		
		if @cycle > 0
			betaI = ((oldgrad[1] - oldgrad[0]) - @h0 * (oldvec[1] - oldvec[0])).abs / (oldvec[1] - oldvec[0]).abs
		else
			betaI = 1.0/@settings[:opt_diagonal_h0]
		end
		
		h0_csr_betaI = SparseMatrixCSX.from_sparse_matrix(:row, @h0)
		
		i = 0
		# repeats until angle direction,gradient is < 90
		begin
			# converts sparse matrix to standard CSR format (due to UMFPack)
			
			(3 * @n_atoms).times{|i|
				h0_csr_betaI[i,i] += betaI
			}
			betaI *= 1.602
			
			start = Time.now
			if h0_csr_betaI.respond_to?(:solve)
				# solves sparse linear system with UMFPack
				z = h0_csr_betaI.solve(q)
			else
				Cuby::raise("To use PLBFGS optimizer, the Algebra extension must be compiled with support for the UMFPACK library")
			end
			stop = Time.now
			Cuby::log.puts_debug "PLBFGS sparse solver timing: #{stop - start} s"
			@stat_counters[":indef_matrix"] if i > 0
			
			i += 1
			Cuby::log.puts_debug "Making matrix PD exceeded 10 steps!!! Giving up." if i > 10
			break if i > 10
		end while Math::acos(z.dot(@gradient) / (z.abs * @gradient.abs)) > Math::PI / 2
		
		i = m - 1
		while i > 0
			beta = ro[i]* y[i].dot(z)
			z += s[i] * (alpha[i] - beta)

			i -= 1
		end
		
		Cuby::log.puts_debug "Angle btwn gradient and direction = #{Math::acos(z.dot(@gradient) / (z.abs * @gradient.abs)) * 180 / Math::PI}"
		z *= -1.0
		
		return z
	end
	
	def step
		# Save current data
		@energy_p = @energy

		# Quque calculation of the point and displacements used for numerical hessian update
		@calculation_queue_proc.call(0,@vector)
		direction_queue

		# This triggers execution of the queue:
		@energy, @gradient = @calculation_run_proc.call(0)
		@stat_counters[:energy_gradient_calls] += 1

		# Energy change
		@delta_e = @energy - @energy_p if @energy_p

		# Print step
		@print_step_proc.call

		# Check convergence
		if @cycle > 0
			return :converged if @check_convergence_proc.call(@delta_e) 
			@stat_counters[:steps_upward] += 1 if @delta_e > 0.0
		end

		# Save coordinates
		@oldgrad.unshift(@gradient)
		@oldvec.unshift(@vector)
		@oldgrad.slice!(@lbfgs_n + 1) # Keep only last N records + 1 so that step can be discarded
		@oldvec.slice!(@lbfgs_n + 1)
		
		# Propose step, update vector
		direction = direction(@oldgrad, @oldvec)
		@step = linesearch(direction)
		
		@vector += @step

		# Success
		return :continue
	end
	
end
