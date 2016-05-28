# Atomic thermochemical data, template for ya YAML data file

class AtomThermochemistry
	attr_reader :h0
	attr_reader :h298_corr
	attr_reader :s0_298

	def initialize(h0, h298_corr, s0_298)
		@h0 = h0
		@h298_corr = h298_corr
		@s0_298 = s0_298
	end
end
