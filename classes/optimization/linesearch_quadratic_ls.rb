module LinesearchQuadraticLS

	def eval(x, direction)
		@calculation_queue_proc.call(0, @vector + x * direction)
		
		# This triggers execution of the queue:
		e, g = @calculation_run_proc.call(0)
		@stat_counters[:energy_gradient_calls] += 1
		return e, g
	end
	
	def linesearch(direction)
		direction = scale_step(direction, @maxstep)
		e1, g1 = eval(1, direction)
		
		if e1 < @energy   # here cannot be "+ 0.1 * @gradient.dot(direction)" so that following conditions were valid
			@reuse_point = [e1, g1]
			@stat_counters[:LS_qls_nth] += 1
			Cuby::log.puts_debug "Linesearch: energy decrease, accept step."
			return direction
		end
		
		d0 = @gradient.dot(direction / direction.abs)
		d1 = g1.dot(direction / direction.abs)
		
		if d0 >= 0   # positive derivative! caused by indefinite preconditioner matrix. doing small step in direction btwn grad and direction (=> downhill)
			@stat_counters[:LS_qls_positive_derivative] += 1
			Cuby::log.puts_debug "Linesearch: derivative in step direction is positive! Doing linear combination of grad and direction."
			dir2 = 0.5 * (direction + @gradient / @gradient.abs * direction.abs)
			return dir2
		end
		
		min = 1
		alpha = 0.8   # váha gradientu ve vých. bodě, 1-alpha v návrhu
		min_part = 0.1
		e0 = @energy
		
		# computing argmin of quadratic function through e0, e1, having least weighted squares from real derivatives (weights are alpha, 1-alpha)
		a = alpha*(d0 - e1 + e0) - (1-alpha)*(d1 - e1 + e0)
		a != 0 ? min = (+e1 - e0 + a) / (2 * a) : min = min_part
		
		if 0.5 < min   # derivative d1 is extremely steep down, quadr. function is concave
			@stat_counters[:LS_qls_extr_steep] += 1
			Cuby::log.puts_debug "Linesearch: very steep descent in proposed point."
			return min_part * direction
		elsif min < 0   # derivative d1 is high, not too much
			@stat_counters[:LS_qls_steep] += 1
			Cuby::log.puts_debug "Linesearch: steep descent in proposed point."
			return min_part * direction
		elsif 0 < min && min < min_part   # similar to the previous case
			@stat_counters[:LS_qls_short_step] += 1
			Cuby::log.puts_debug "Linesearch: proposed step shorter than #{min_part} direction."
			return min_part * direction
		else
			@stat_counters[:LS_qls_ok] += 1
			Cuby::log.puts_debug "Linesearch: OK, uses #{min} direction"
			return scale_step(min * direction, @maxstep)
		end
	end
end
