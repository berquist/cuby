class CpeParameters
	attr_reader :z, :b, :rl, :ru

	def initialize(z, b, rl, ru)
		@z, @b, @rl, @ru = z, b, rl, ru
	end

	def to_hsd(element, indent)
		s =  " " * indent + element.to_s + " {\n"
		s += " " * indent + "  Z = #{@z}\n"
		s += " " * indent + "  B = #{@b}\n"
		s += " " * indent + "  Rl = #{@rl}\n"
		s += " " * indent + "  Ru = #{@ru}\n"
		s += " " * indent + "}\n"
	end

end
