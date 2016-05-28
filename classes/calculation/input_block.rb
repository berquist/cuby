# Container for the description of input blocks

class InputBlock
	attr_accessor :name
	attr_accessor :status
	attr_accessor :description

	def initialize(name, status, description)
		@name = name
		@status = status
		@description = description
	end

	def InputBlock.[](name, status, description)
		return InputBlock.new(name, status, description)
	end
end
