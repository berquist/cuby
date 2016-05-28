################################################################################
#
# Class ExpressionParser
#
# Author: Jan Rezac
# Date created: 2010-04-16
# License: Cuby license
# Description: Parser of infix notation expressions
# Status: Tested, Documented
#
################################################################################

class ExpressionParser

	def initialize
		@operators = []
		@functions = []
		@actions = {}

		@parenthesis = ["(",")"]
	end

	def leaf_processor(action_proc)
		#: A Proc that takes the leaf of the tree (the value) and converts it to
		#: the form the parser works on (i.e. convert string to number)
		@actions[:value] = action_proc
	end

	def add_operator(character, priority, action_proc)
		#: Define new operator by its name (one character), priority level
		#: and a Proc that executes it on left and right value

		# Add priority level if it was missing
		unless @operators[priority]
			@operators[priority] = []
		end

		# Add operator to operators list
		@operators[priority] << character

		# Add action
		@actions[character] = action_proc
	end

	def add_function(name, action_proc)
		#: Define new function of a given name, the parameters are pased to
		#: the Proc as a string.
		@functions << name
		@actions[name] = action_proc
	end

	def parse(string)
		#: Parse a string, return the value

		# Look for operators, high priority first
		@operators.each{|operator_set| # Iterate through operator priorities
			depth = 0
			string.size.times{|i|
				c = string[i..i] # Ruby 1.8 compatible way to get a character
				if c == @parenthesis[0]
					depth += 1
					next
				end
				if c == @parenthesis[1]
					depth -= 1
					next
				end
				if depth == 0 && operator_set.include?(c)
					l = string[0..i-1]
					r = string[i+1..-1]
					raise("Expression can not start with an operator") if l == ""
					raise("Expression can not end with an operator") if r == ""
					return @actions[c].call(parse(l),parse(r))
				end
			}
			if depth != 0
				raise("Parenthesis not matching in expression")
			end
		}
		# Parenthesis
		if string[0..0] == @parenthesis[0] && string[-1..-1] == @parenthesis[1]
			return parse(string[1..-2])
		end
		
		# Functions
		@functions.each{|func|
			if match = /^#{func}\((.*)\)$/.match(string)
				return @actions[func].call(match[1])
			end
		}

		# Leaf: a value
		return @actions[:value].call(string)
	end

end

__END__
# Tests
parser = ExpressionParser.new
parser.leaf_processor(Proc.new{|str| str.split(",").map{|x| x.to_i}})
parser.add_operator("|", 1, Proc.new{|l,r| l | r})
parser.add_operator("&", 0, Proc.new{|l,r| l & r})
parser.add_function("%seq", Proc.new{|parameters| a,b = parameters.split(";"); (a.to_i..b.to_i).to_a})

puts parser.parse(ARGV[0])

