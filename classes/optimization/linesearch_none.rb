module LinesearchNone
	
	def linesearch(direction)
		return scale_step(direction, @maxstep)
	end
	
end