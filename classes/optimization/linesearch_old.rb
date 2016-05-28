module LinesearchOld

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
		
		if e1 < @energy + 0.1 * @gradient.dot(direction)   # no guaranties, can be here
			@reuse_point = [e1, g1]
			@stat_counters[:LS_nth] += 1
			Cuby::log.puts_debug "Linesearch: enough energy decrease, accept step."
			return direction
		end
		
		d0 = @gradient.dot(direction / direction.abs)
		d1 = g1.dot(direction / direction.abs)
		
		min = 1
		min = - d0 / (d1 - d0) if d1 != d0
		
		if 0.618 < min && min < 1.618
			@reuse_point = [e1, g1]
			@stat_counters[:linesearch_btwn_0_618__1_618] += 1
			return direction
		elsif 0 < min && min < 0.382
			@stat_counters[:linesearch_less_0_382] += 1
			return 0.382 * direction
		elsif min < 0
			@stat_counters[:linesearch_negative] += 1
			return 0.382 * direction
		elsif min > 2.618
			@stat_counters[:linesearch_more_2_618] += 1
			return scale_step(2.618 * direction, @maxstep)
		else
			@stat_counters[:linesearch_ok] += 1
			return scale_step(min * direction, @maxstep)
		end
	end
end