module LinesearchCubic

	def eval(x, direction)
		@calculation_queue_proc.call(0, @vector + x * direction)
		
		# This triggers execution of the queue:
		e, g = @calculation_run_proc.call(0)
		
		@stat_counters[:energy_gradient_calls] += 1
		return e, g
	end
	
	def linesearch(direction)
		direction = scale_step(direction, @maxstep)
		
		e0 = @energy
		e1, g1 = eval(1, direction)
		d0 = @gradient.dot(direction)   # function scaled to interval [0,1], derivative must be multiplied by direction.abs
		d1 = g1.dot(direction)
		
		a =  2*e0 +   d0 - 2*e1 + d1
		b = -3*e0 - 2*d0 + 3*e1 - d1
		c =           d0
		d =    e0
		
		#~ Cuby::log.puts_debug "#{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"   # interpolation polynomial
		#~ Cuby::log.puts_debug "e0=#{e0}, d0=#{d0}, e1=#{e1}, d1=#{d1}"
		
		#~ #!# smazat
		#~ if @cycle == 2 or @cycle == 4 or @cycle == 11 or @cycle == 13
			#~ n = 20
			#~ n.times do |i|
				#~ ee, gg = eval((i+1).to_f/n, direction)
				#~ puts "Energ#{@cycle}: #{(i+1).to_f/n}\t#{ee}"
			#~ end
		#~ end
		
		
		min_part = 0.1
		max_part = 2.0
		min_effic = 2.0
		max_trust_der = 10.0
		
		if d0 >= 0   # positive derivative! possibly caused by indefinite preconditioner matrix. doing small step in direction btwn grad and direction (=> downhill)
			@stat_counters[:LS_positive_derivative] += 1
			Cuby::log.puts_debug "Linesearch: Positive derivative! The angle is #{Math::acos(-@gradient.dot(direction)/(@gradient.abs * direction.abs)) * 180/Math::PI} deg."
			#!# rather use some history flush method !!!
			@oldgrad = [@oldgrad[0],@oldgrad[1]]
			@oldvec = [@oldvec[0],@oldvec[1]]
			Cuby::log.puts_debug "#{@energy} Energy:"
			return min_part * (direction + @gradient / @gradient.abs * direction.abs)
		end
		
		if (b**2 - 3*a*c) < 0   # minimum does not exist
			@stat_counters[:LS_REUSE_no_minimum] += 1
			Cuby::log.puts_debug "Linesearch: No minimum. #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
			@reuse_point = [e1, g1]
			return direction
		else   # minimum does exist
			if d1.abs / d0.abs > max_trust_der   # new derivative is extremely high, doing smallest possible step
				@stat_counters[:LS_extreme_derivative] += 1
				Cuby::log.puts_debug "Linesearch: Extreme derivative! d0 = #{d0}, d1 = #{d1}; #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
				Cuby::log.puts_debug "#{@energy} Energy:"
				return min_part * direction
			end
			
			mins = [(-b - (b**2 - 3*a*c)**0.5) / (3*a), (-b + (b**2 - 3*a*c)**0.5) / (3*a)].sort   # [].sort because 'a' can be negative
			mins[0] > 0 ? min = mins[0] : min = mins[1]   # min = smaller positive value from mins, must be in (0,1)
			
			if min <= 0
				@stat_counters[:LS_REUSE_min_less_0] += 1
				Cuby::log.puts_debug "Linesearch: Min < 0. #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
				@reuse_point = [e1, g1]
				return direction
			elsif 0 < min and min < 1
				if min < min_part
					@stat_counters[:LS_min_less_min_part] += 1
					Cuby::log.puts_debug "Linesearch: Min near 0. #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
					Cuby::log.puts_debug "#{@energy} Energy:"
					return min_part * direction
				elsif min > 1 - min_part
					@stat_counters[:LS_REUSE_min_near_1] += 1
					Cuby::log.puts_debug "Linesearch: Min near 1. #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
					@reuse_point = [e1, g1]
					return direction
				elsif e1 < e0 and min_effic*(e0-e1).abs > (e0-(a*min**3+b*min**2+c*min+d)).abs
					@stat_counters[:LS_REUSE_small_estimated_benefit] += 1
					Cuby::log.puts_debug "Linesearch: Small benefit (0,1). #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
					@reuse_point = [e1, g1]
					return direction
				elsif e1 >= e0
					@stat_counters[:LS_energy_increase] += 1
					Cuby::log.puts_debug "Linesearch: Energy increase. #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
					Cuby::log.puts_debug "#{@energy} Energy:"
					return min * direction
				else
					@stat_counters[:LS_good_benefit] += 1
					Cuby::log.puts_debug "Linesearch: Good benefit (0,1). #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
					Cuby::log.puts_debug "#{@energy} Energy:"
					return min * direction
				end
			else   # min >= 1
				if min < 1 + min_part
					@stat_counters[:LS_REUSE_min_near_1] += 1
					Cuby::log.puts_debug "Linesearch: Min near 1. #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
					@reuse_point = [e1, g1]
					return direction
				else   # min >= 1+min_part
					min = [min, max_part].min
					if min_effic*(e0-e1).abs > (e0-(a*min**3+b*min**2+c*min+d)).abs
						@stat_counters[:LS_longer_step] += 1
						Cuby::log.puts_debug "Linesearch: Good benefit (1, ...). #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
						Cuby::log.puts_debug "#{@energy} Energy:"
						return min * direction   #!# use scale_step ???
					else
						@stat_counters[:LS_REUSE_small_estimated_benefit] += 1
						Cuby::log.puts_debug "Linesearch: Small benefit (1, ...). #{a}*x**3 + #{b}*x**2 + #{c}*x + #{d}"
						@reuse_point = [e1, g1]
						return direction
					end
				end
			end
			
		end
		
	end
end
