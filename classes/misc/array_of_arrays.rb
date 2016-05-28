class ArrayOfArrays < Array
	# An array of arrays of different lengths

	def no_of_combinations
		# Return the number of all possible combinations of the items in the iner arrays
		total = 1
		each{|a|
			total *= a.size
		}
		return total
	end

	def iterate_all_combinations
		# Iterate over all possible combinations of the items in the iner arrays
		# Yields an array of indexes
		no_of_combinations.times{|i|
			indexes = []
			x = i
			self.size.times{|j|
				divider = self[self.size - 1 - j].size
				indexes.unshift(x % divider)
				x = x / divider
			}
			yield indexes
		}

		return nil
	end
end

