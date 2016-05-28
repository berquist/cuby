
module OptimizerCommon

	def tr_update
		#: Trust radius update

		return :disabled if @settings[:opt_tr_update] == false

		# Update is not applied in cycle 0
		return :skipped if @cycle < 1

		r = @delta_e / @predicted_de
		if r > 0.75
			@maxstep *= 2.0
			status = :step_increased
		end
		if r < 0.25
			if @maxstep > @step.abs
				# If maxstep was not applied, but step was refused, set maxstep to current step size
				@maxstep = @step.abs
			end
			@maxstep /= 3.0
			status = :step_decreased
		end
		if @maxstep > @maxstep_max
			@maxstep = @maxstep_max 
			status = :step_maximum
		end
		if @maxstep < @maxstep_min
			@maxstep = @maxstep_min 
			status = :step_minimum
		end

		if r < 0 && @settings[:opt_refuse_steps]
			# Restore previous step
			@energy = @energy_p
			@vector = @vector_p
			@gradient = @grad_p
			@hessian = @hessian_p

			status = :refused
			@stat_counters[:refused_steps] += 1

			# Failure - defined as second bad step with step size limited
			# to minimum
			if @maxstep == @maxstep_min
				@failure_count += 1
				if  @failure_count >= 2
					#!# Optionally, reset hessian on first failure and stop the optimization
					#!# only on a second one.
					status = :failed
				end
			end
		else
			@failure_count = 0
		end
		Cuby.log.puts_debug "Step: #{status} (r = #{r})"
		return status
	end

	def scale_step(s, maxstep)
		scaled = false
		case @step_scaling
		when :abs
			if s.abs > maxstep
				 s = s / s.abs * maxstep
				 scaled = true
			end
		when :max
			if s.max > maxstep
				 s = s / s.max * maxstep
				 scaled = true
			end
		when :components
			s.each_index{|i|
				if s[i].abs > maxstep
					s[i] = maxstep * s[i]/s[i].abs
					scaled = true
				end
			}
		else
			raise "Unknown step scaling mode '#{@step_scaling}'"
		end

		if scaled
			 Cuby::log.puts_debug "Maxstep: Step size limited"
		else
			 Cuby::log.puts_debug "Maxstep: Step size OK"
		end

		return s
	end

	def hessian_update_qn
		# Update is not applied in cycle 0
		return if @cycle < 1

		dx_v = @step
		y_v = (@gradient - @grad_p)

		dx = dx_v.to_matrix
		y = y_v.to_matrix

		y_dot_dx =  y_v * dx_v


		# Skip hessian update in bad cases to conserve positive definite hessian, 
		# also accounting for numerical precision (Dennis & Schnabel 1984)
		unless y_dot_dx > Math::sqrt(1.0e-10) * dx_v.abs * y_v.abs
			return
		end

		#--------------------------------------------------
		# if cycle == 1 && @update_diagonal_hessian
		#	 # Update diagonal estimate
		#	 diagonal = y_v.dot(dx_v) / y_v.abs**2
		#	 @hessian.n.times{|i|
		#		 @hessian[i,i] = 1.0/diagonal
		#	 }
		# end
		#-------------------------------------------------- 

		case @settings[:opt_qn_update_formula]
		when :bfgs # Broyden-Fletcher-Shanno-Goldfarb
			@hessian = @hessian + 
				(y * y.transpose)/y_dot_dx -
				(@hessian * dx)*(@hessian*dx).transpose /
				(dx.transpose * @hessian * dx)[0,0]
		when :dfp # Davidon-Fletcher-Powell
			i = Matrix.unit(@hessian.m)
			@hessian = (i - y*dx.transpose / y_dot_dx) * 
				@hessian * 
				(i - dx*y.transpose / y_dot_dx) +
				(y * y.transpose)/y_dot_dx
		when :broyden
			@hessian = @hessian +
				(y - @hessian*dx)/(dx.transpose * dx)[0,0] *
				dx.transpose
		when :sr1
			@hessian = @hessian +
				((y - @hessian*dx) * (y - @hessian*dx).transpose) /
				((y - @hessian*dx).transpose * dx)[0,0]
		when :dfp_bfgs
			y_dot_dx = y_dot_dx
			### The calculation of S is taken from the formula for inverse hessian, there might be a more efficient 
			### way to do it (withou the inverse)
			s = 1.0 + (y.transpose * @hessian.inverse * y)[0,0] / y_dot_dx
			if 1.0 < (s - 1.0) && y_dot_dx > 0.0
				# DFP
				i = Matrix.unit(@hessian.m)
				@hessian = (i - y*dx.transpose / y_dot_dx) * 
					@hessian * 
					(i - dx*y.transpose / y_dot_dx) +
					(y * y.transpose)/y_dot_dx
			else
				# BFGS
				@hessian = @hessian + 
					(y * y.transpose)/y_dot_dx -
					(@hessian * dx)*(@hessian*dx).transpose /
					(dx.transpose * @hessian * dx)[0,0]
			end
		when :sr1_bfgs
			# SR1 and BFGS mixing using Bofill factor
			# This is used in Gaussian - DOI: 10.1021/ct050275a

			# SR1 update matrix
			sr1_upd = ((y - @hessian*dx) * (y - @hessian*dx).transpose) /
				((y - @hessian*dx).transpose * dx)[0,0]

			# BFGS update matrix
			bfgs_upd = (y * y.transpose)/y_dot_dx -
				(@hessian * dx)*(@hessian*dx).transpose /
				(dx.transpose * @hessian * dx)[0,0]

			# Bofill factor
			hdx_y = @hessian*dx - y
			factor = ( (hdx_y.transpose * dx)[0,0]**2 /
				((hdx_y.transpose * hdx_y)[0,0] * (dx.transpose*dx)[0,0]))**0.5

			# Perform the update
			@hessian = @hessian + factor * sr1_upd + (1.0-factor) * bfgs_upd
		else
			raise("Invalid hessian update formula")
		end

	end

	def hessian_update_inverse_qn
		# Update is not applied in cycle 0
		return if @cycle < 1

		dx_v = @step
		y_v = (@gradient - @grad_p)

		dx = dx_v.to_matrix
		y = y_v.to_matrix

		y_dot_dx =  y_v * dx_v

		# Skip hessian update in bad cases to conserve positive definite hessian, 
		# also accounting for numerical precision (Dennis & Schnabel 1984)
		unless y_dot_dx > Math::sqrt(1.0e-10) * dx_v.abs * y_v.abs
			return
		end

		# Cycle 1: diagonal update 
		#!# check whether it works
		#--------------------------------------------------
		# if @cycle == 1
		#	 # Update diagonal estimate
		#	 diagonal = y_v.dot(dx_v) / y_v.abs**2
		#	 @hessian.n.times{|i|
		#		 @hessian[i,i] = diagonal
		#	 }
		# end
		#-------------------------------------------------- 

		case @settings[:opt_qn_update_formula]
		when :bfgs # Broyden-Fletcher-Shanno-Goldfarb
			a = Matrix.unit(@hessian.m) - (y * dx.transpose) / y_dot_dx
			@hessian = a.transpose * @hessian * a + (dx * dx.transpose) / y_dot_dx
		when :dfp # Davidon-Fletcher-Powell
			@hessian += (dx * dx.transpose) / y_dot_dx -
				@hessian * y * y.transpose * @hessian.transpose / (y.transpose * @hessian * y)[0,0]
		when :broyden
			@hessian += (dx - @hessian * y) * y.transpose * @hessian / (y.transpose * @hessian * y)[0,0]
		when :sr1
			a = (dx - @hessian * y)
			@hessian += a * a.transpose / (a.transpose * y)[0,0]
		when :dfp_bfgs
			y_dot_dx = y_dot_dx
			s = 1.0 + (y.transpose * @hessian * y)[0,0] / y_dot_dx
			if 1.0 < (s - 1.0) && y_dot_dx > 0.0
				# DFP
				@hessian += (dx * dx.transpose) / y_dot_dx -
					@hessian * y * y.transpose * @hessian.transpose / (y.transpose * @hessian * y)[0,0]
			else
				# BFGS
				a = Matrix.unit(@hessian.m) - (y * dx.transpose) / y_dot_dx
				@hessian = a.transpose * @hessian * a + (dx * dx.transpose) / y_dot_dx
			end
		when :sr1_bfgs
			# SR1 and BFGS mixing using Bofill factor
			# This is used in Gaussian - DOI: 10.1021/ct050275a

			# SR1 hessian
			a = (dx - @hessian * y)
			sr1_hessian = @hessian + a * a.transpose / (a.transpose * y)[0,0]

			# BFGS hessian
			a = Matrix.unit(@hessian.m) - (y * dx.transpose) / y_dot_dx
			bfgs_hessian = a.transpose * @hessian * a + (dx * dx.transpose) / y_dot_dx

			# Bofill factor
			#!# Formula taken from non-inverse hessian form requires matrix inverse
			hdx_y = @hessian.inverse*dx - y
			factor = ( (hdx_y.transpose * dx)[0,0]**2 /
				((hdx_y.transpose * hdx_y)[0,0] * (dx.transpose*dx)[0,0]))**0.5

			# Perform the update
			@hessian = factor * sr1_hessian + (1.0-factor) * bfgs_hessian
		else
			raise "Invalid hessian update formula"
		end
	end
end
