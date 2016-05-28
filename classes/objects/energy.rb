# Notes on Energy class design:
# The construction of classes for different types of energy and their
# components is complicated, but achieves two important goals:
# 1) Type of the energy is printed along the value
# 2) Type of the energy is printed along the decomposition
# 3) Nice and short YAML form

# TBD?
# * support for short yaml notation when no components are present, something like 'object/energy: 1.00'
# * Global parameter for number format read from Cuby::Log
# * Verbosity support (automated printing based on verbosity), e.g. method Energy.print_log
#

#===============================================================================
# Various kinds of Energy and EnergyComponents are generated from the prefixes
# (the Energy class body comes later)
#===============================================================================

class Energy
	PREFIXES = [
		'Interaction'
	]
end

class EnergyComponents < Hash
end

# Class generator for energies with various prefixes
Energy::PREFIXES.each{|prefix|
	eval "class #{prefix}Energy < Energy; end"
	eval "class #{prefix}EnergyComponents < EnergyComponents; end"
}

#===============================================================================
# Energy class
#===============================================================================

class Energy
	attr_reader :value
	attr_reader :components

	def initialize(value, components = nil)
		@value = value
		set_components(components) unless components.nil?
	end


	def to_s(options = {})
		opts = {
			:format => "%.4f",
			:prefix => "#{energy_type_name}: ",
			:unit => true
		}
		opts.merge!(options)
		s = opts[:prefix] + sprintf(opts[:format],@value)
		s += ' kcal/mol' if opts[:unit]
		return s
	end

	def set_components(hash)
		@components = eval("#{self.class.to_s}Components").new(hash)
	end

	def energy_type_name
		return self.class.to_s.gsub(/([a-z])([A-Z])/,'\1 \2').downcase.capitalize
	end

end

#===============================================================================
# EnergyComponents class
#===============================================================================

class EnergyComponents < Hash

	ALLOWED_NAMES = [
		# HF calculations
		:scf_energy,

		# Correlated calculations
		:mp2_correlation,
		:mp3_correlation,
		:mp4_correlation,
		:ccsd_correlation,
		:"ccsd[t]_correlation",
		:"ccsd(t)_correlation",
		:ccsdt_correlation,
		:"ccsdt[q]_correlation",
		:"ccsdt(q)_correlation",
		:ccsdtq_correlation,

		# Modifiers
		:dispersion_correction,
		:hbond_correction,

		# Forcefield
		:bonds,
		:angles,
		:torsions,
		:lennard_jones,
		:electrostatic
	]

	def initialize(hash = {})
		hash.each_pair{|k,v| self[k] = v}
	end

	def []=(index, value)
		check_key(index)
		super(index, value)
	end

	def to_s(options = {})
		opts = {
			:format => "%-25s%.4f",
			:prefix => "#{energy_type_name} decomposition: ",
			:unit => true,
			:separator => "\n",
			:indent => "   "
		}
		opts.merge!(options)

		unit = ''
		unit = ' kcal/mol' if opts[:unit]

		a = []
		a << opts[:prefix]
		each_pair{|key, value|
			a << opts[:indent] + sprintf(opts[:format], key.to_s + ':', value) + unit
		}

		return a.join(opts[:separator])
	end

	def energy_type_name
		return self.class.to_s.gsub(/Components$/,'').gsub(/([a-z])([A-Z])/,'\1 \2').downcase.capitalize
	end

	def check_key(key)
		raise "Energy component name must be a Symbol" unless key.class == Symbol
		raise "Energy component name '#{key}' is not on the list of alloved names" unless ALLOWED_NAMES.include?(key)
		return true
	end
end

#===============================================================================
# Testing
#===============================================================================

require "yaml"

e = Energy.new(1.0)
puts e


e = InteractionEnergy.new(1.0, {:scf_energy => 0.5, :mp2_correlation => 0.5})
puts e
puts e.components
puts

puts [e].to_yaml
