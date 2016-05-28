#===============================================================================
# Lookup table of binomial coefficients
#===============================================================================

module BinomialCoefficients
	# Initialize with max. N = 100
	DEFAULT_MAX = 100

	# Module data
	@@f_by_n = []
	@@max = 0

	#=======================================================================
	# Access to the coefficients
	#=======================================================================

	def BinomialCoefficients.f(k,n)
		# Binomial coefficient
		# in notation (n over k) or F_k(n)
		raise "Binomial coefficient n is outside of the lookup table" if n > @@max
		raise "Binomial coefficent outside of range, k,n = #{k}, #{n}" if @@f_by_n[n][k] == nil
		return @@f_by_n[n][k]
	end

	def BinomialCoefficients.f0(k,n)
		# Binomial coefficent or zero if k > n
		return 0 if k > n
		return BinomialCoefficients.f(k,n)
	end

	def BinomialCoefficients.nm(n,m)
		# Binomial coefficient
		raise "Binomial coefficient n is outside of the lookup table" if n > @@max
		return 0 if m > n
		raise "Binomial coefficent outside of range, m,n = #{m}, #{n}" if @@f_by_n[n][m] == nil
		return @@f_by_n[n][m]
	end

	#=======================================================================
	# Generators
	#=======================================================================
	
	def BinomialCoefficients.generate_f(max)
		# Generate binomial coefficients
		
		@@max = max
	
		# Recursive formula avoiding factorials	
		@@f_by_n = []
		(0..max).each{|n|
			@@f_by_n[n] = []
			@@f_by_n[n][0] = 1
			@@f_by_n[n][n] = 1
			(1..n-1).each{|k|
				@@f_by_n[n][k] = @@f_by_n[n-1][k-1] + @@f_by_n[n-1][k]
			}
		}

		# Test: prints Pascal's triangle
		#@f_by_n.each{|a_n| puts a_n.join(" ")}
		return nil
	end 

	# Module initialization
	BinomialCoefficients.generate_f(DEFAULT_MAX)
end

