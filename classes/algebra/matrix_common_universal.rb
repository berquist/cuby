################################################################################
#
#  module MatrixCommonUniversal : Common methods for Matrix and SparseMatrix
#
################################################################################

module Algebra
module MatrixCommonUniversal

	# Prerequisites to include this module:
	# .m, .n
	# .row_as_array(i)

	#=======================================================================
	# basics
	#=======================================================================
	
	def dimensions
		#: Dimensions of the matrix as [m,n]
		return [m,n]
	end

	def size
		#: Dimensions of the matrix as [m,n]
		return [m,n]
	end

	def size_str
		# Dimensions of the matrix as a string "m x n"
		return "#{m} x #{n}"
	end

	#=======================================================================
	# type conversion
	#=======================================================================
	
	def to_s(format="%10.6f", separator = " ", decoration_s = "| ", decoration_e = " |")
		#: Converts matrix to string. Format parameter can be used to modify
		#: the numbers.
		s = ''
		rows = []

		self.m.times {|i|
			line = row_as_array(i).collect{|f| 
				f = 0.0 if f == -0.0
				sprintf(format,f)
			}
			rows << decoration_s + line.join(separator) + decoration_e
		}

		return rows.join("\n")
	end

	def inspect(format="%5.2f")
		#: Equal to -->to_s with prepended newline for nice irb output.
		#: Default format is shorter for convenience.
		return "\n" + to_s(format)
	end

	def to_s_rb(format="%5.2f")
		#: Converts matrix to string. Format parameter can be used to modify
		#: the numbers.
		s = ''
		rows = []

		self.m.times {|i|
			line = row_as_array(i).collect{|f| 
				f = 0.0 if f == -0.0
				sprintf(format,f)
			}
			rows << "       [#{line.join(", ")}]"
		}

		rows.first.gsub!("       ","Matrix[")
		rows.last << "]"
		s = rows.join(",\n")
	end

	def to_s_triangular(format="%10.6f", separator = " ")
		#: Converts matrix to string, printing only half (j<=i) of the matrix. Format parameter can be used to modify
		#: the numbers.
		s = ''
		rows = []

		self.m.times {|i|
			line = []
			(i + 1).times{|j|
				line << self[i,j]
			}
			line.map!{|f| 
				f = 0.0 if f == -0.0
				sprintf(format,f)
			}
			rows << line.join(separator)
		}

		return rows.join("\n")
	end

end
end
