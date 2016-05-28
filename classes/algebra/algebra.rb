################################################################################
#
# Linear Algebra Library
#
# Author: Jan Rezac
# Date created: 2008-11-03
# License: Cuby license
# Description: Header file for algebra classes Matrix and Vector
# Status: Tested, Documented
#
################################################################################

#: This is the module to be included for use of Vector and Matrix classes.
module Algebra
	@@directory = File.dirname(__FILE__)

	# Load classes
	if (FileTest.exists?(@@directory + "/algebra_c.so") || FileTest.exists?(@@directory + "/algebra_c.bundle")) && !$algebra_force_ruby
		# Ruby version match test
		compiled_version = IO.readlines(File.dirname(__FILE__) + "/_version")[0].strip
		if compiled_version != RUBY_VERSION
			raise "Extension 'algebra' had been compiled with another version of ruby"
		end
		# Use binary if possible
		require @@directory+"/algebra_c"
		require @@directory+"/matrix_c_fixes"
		# When compiled without lapack, load equivalent ruby routines
		require @@directory+"/matrix_ruby_lapack" unless Matrix.method_defined?(:eigensystem)
	else
		# Otherwise, load ruby equivalent
		require @@directory+"/matrix_ruby"
		require @@directory+"/matrix_ruby_lapack"
		require @@directory+"/vector_ruby"
		require @@directory+"/sparse_matrix_csx/sparse_matrix_csx_ruby"
		# Make the module respond to version query
		def version; return :ruby; end
	end

	# Common parts, in ruby
	require @@directory+"/matrix_common"
	require @@directory+"/vector_common"
	require @@directory+"/sparse_matrix_csx/sparse_matrix_csx_common"

	# Sparse matrix, in ruby (no C equivalent)
	require @@directory+"/sparse_matrix/sparse_matrix.rb"
end




################################################################################

#===============================================================================
# Numeric classes extension
#===============================================================================

#: Multiplication operators of Float and Fixnum are extended to handle any
#: multiplication of objects that respond to method "comutative_multiply"
#: by calling Something * Number when the actual code is Number * Something.

class Float
	alias mul_without_extension *

	def *(operand) # => usually instance of operand.class
		#: Multiplication is extended to handle Float * operand
		if operand.respond_to?(:comutative_multiply)
			return operand.comutative_multiply(self)
		else
			return mul_without_extension(operand)
		end
	end

end

class Fixnum
	alias mul_without_extension *

	def *(operand) # => usually instance of operand.class
		#: Multiplication is extended to handle Fixnum * operand
		if operand.respond_to?(:comutative_multiply)
			return operand.comutative_multiply(self)
		else
			return self.mul_without_extension(operand)
		end
	end

end
