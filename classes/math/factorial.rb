#===============================================================================
# Fast factorial evaluation using a lookup table
#===============================================================================

module Factorial
	# Initialize with max. N = 100
	# Can be extended later by calling Factorial.add(N)
	DEFAULT_MAX = 100

	# Module data
	@@factorials = []
	@@max = 0

	def Factorial.n(x)
		# Returns the factorial
		if x > @@max
			raise "Range of pre-generated factorials exceeded with n = #{x}"
			# Here, the number of stored factorials can be increased on the go
			# by calling the add method, but some sensible limits should be checked
		end
		return @@factorials[x]
	end

	def Factorial.add(max)
		# This adds factorials up to max to the storage

		# If the array is empty, initialize the first one
		if @@factorials.size == 0
			@@factorials[0] = 1
		end
		# Start the cycle with the last one present
		x = @@factorials.size
		f = @@factorials.last
		# Add more
		while x <= max
			f *= x
			@@factorials[x] = f
			x += 1
		end
		# Save size
		@@max = @@factorials.size
		return nil
	end

	# Module initialization
	Factorial.add(DEFAULT_MAX)
end

