class DotWriter
	def initialize(line_length = 80)
		@line_length = line_length
		@counter = 0
	end

	def print_dot(dot=".")
		print(dot)
		@counter += 1
		if @counter == @line_length
			@counter = 0
			print("\n")
		end
		return nil
	end

	def finish
		print("\n") unless @counter == 0
	end
end
