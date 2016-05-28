class ChemicalReaction
	attr_reader :reactants
	attr_reader :products

	class CRItem
		attr_reader :name
		attr_reader :num

		attr_accessor :job # Link to the calculation

		def initialize(num, name)
			@num = num
			@name = name
		end

		def to_s
			return "#{@num} #{@name}"
		end
	end

	def initialize(reactants, products)
		@reactants = reactants
		@products = products
	end

	def ChemicalReaction.from_string(s)
		# Get left and right-hand side
		raise "Chemical reaction formula should be separated with '->'" unless s =~ /->/
		left, right = s.split('->')
		# Convert strings into lists of reactants and products
		reaction = ChemicalReaction.new(ChemicalReaction.process_side(left), ChemicalReaction.process_side(right))
		return reaction
	end

	def ChemicalReaction.process_side(s)
		list = []
		s.split('+').each{|item|
			item = item.strip
			next if item == ""
			if item =~ /^[0-9]/
				matchdata = /^([0-9]+)\s*(.*)/.match(item)
				list << CRItem.new(matchdata[1].to_i,matchdata[2])
			else
				list << CRItem.new(1,item)
			end
		}
		return list
	end

	def to_s
		return @reactants.map{|x| x.to_s}.join(" + ") + " -> " + @products.map{|x| x.to_s}.join(" + ")
	end

	def all_items
		return @reactants | @products
	end

end


