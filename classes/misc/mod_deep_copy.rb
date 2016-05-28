################################################################################
#
#  Deep copy
#
#  Duplicates recursively entire object, providing all instance variables could
#  be copied directly, using dup or deep_copy. Array and Hash classes are
#  extended with deep_copy method for their duplication.
#
#  Special treatment is required for some classes that have dup method but it
#  can not be used. Nil, true, false, Fixnum, Bignum, Float and Complex are
#  already fixed, other problematic classes could be added in analogous way.
#
################################################################################

class DeepCopyError < RuntimeError
end

module DeepCopy

	def deep_copy
		# locking support
		if respond_to?(:locked?)
			if self.locked?
				raise "deep_copy can not work on locked objects"
			end
		end

		# deep_copy / dup all instance variables
		result = self.class.new
		ivs = instance_variables
		ivs.each do |iv|
			itv = instance_variable_get(iv)
			otv = case
			      when itv.nil?
				      nil
			      when itv.class == TrueClass
				      true
			      when itv.class == FalseClass
				      false
			      when itv.class == Symbol
				      itv
			      when itv.class == Fixnum || itv.class == Bignum || itv.class == Float #|| itv.class == Complex
				      itv
			      when itv.respond_to?(:deep_copy)
				      itv.deep_copy
			      when itv.respond_to?(:dup)
				      itv.dup
			      else
				      raise(DeepCopyError, "class #{itv.class} can not be copied by deep_copy")
			      end
			result.instance_variable_set(iv, otv)
		end

		return result
	end
end

#===============================================================================
# Manually fix some problematic classes
#===============================================================================

class Array
	def deep_copy
		result = []
		each_index {|i|
			itv = self[i]
			result[i] = case
			      when itv.nil?
				      nil
			      when itv.class == TrueClass
				      true
			      when itv.class == FalseClass
				      false
			      when itv.class == Symbol
				      itv
			      when itv.class == Fixnum || itv.class == Bignum || itv.class == Float #|| itv.class == Complex
				      itv
			      when itv.respond_to?(:deep_copy)
				      itv.deep_copy
			      when itv.respond_to?(:dup)
				      itv.dup
			      else
				      raise(DeepCopyError, "class #{itv.class} can not be copied by deep_copy")
			      end
		}
		return result
	end
end

class Hash
	def deep_copy
		result = {}
		each_key {|key|
			itv = self[key]
			result[key] = case
			      when itv.nil?
				      nil
			      when itv.class == TrueClass
				      true
			      when itv.class == FalseClass
				      false
			      when itv.class == Symbol
				      itv
			      when itv.class == Fixnum || itv.class == Bignum || itv.class == Float #|| itv.class == Complex
				      itv
			      when itv.respond_to?(:deep_copy)
				      itv.deep_copy
			      when itv.respond_to?(:dup)
				      itv.dup
			      else
				      raise(DeepCopyError, "class #{itv.class} can not be copied by deep_copy")
			      end
		}
		return result
	end
end

