################################################################################
#
# Class Atom
#
# Author: Jan Rezac
# Date created: 2008-09-05
# License: Cuby license
# Description: Atom object
# Status: Working, documented
#
################################################################################

# Various properties of the atom are stored in hash properties. There are no
# default value, only the properties already set can be used, and none are
# set by default. It is simply a space to store some temporary information on
# the atom.

require "classes/coordinate/coordinate"
require "classes/misc/mod_deep_copy"
require "classes/periodic_table/periodic_table"

class Atom < Coordinate

	# List of privately used instance variables:
	# @x, @y, @z - inherited from coordinate - should not be used now,
	# C version of coordinate implements these as methods.
	# Use self.x ... instead
	
	# List of public attributes
	attr_accessor 	:element	# Symbol
	attr_reader	:properties	# Properties class(symbol indexed)
	# use set_property to set this variable
	
	ALLOWED_PROPERTIES = {
		# General
		:mass => :float,		# Atomic mass
		:charge => :float,		# Atomic charge
		:dummy => :boolean,		# Boolean flag for dummy 
		:ghost => :boolean,		# Boolean flag for ghost atom (basis set without electrons)

		# Connectivity-derived, defined in geometry_connectivity.rb
		#--------------------------------------------------
		# :hybridization,
		# :bound_atoms,	# Array of indexes of bound atoms
		# :atom_type,
		#-------------------------------------------------- 

		:file_index => :integer, 	# Atom index in the original geometry file (set on file reading, not written on output)

		# Information from PDB file, defined in geometry_readwrite.rb
		:pdb_atom_no => :integer,
		:pdb_atom_name => :string,
		:pdb_res_name => :string,	
		:pdb_res_no => :integer,
		:terafter => :boolean,
		:pdb_heteroatom => :boolean,
		:pdb_chain => :string,
		:pdb_alt_loc => :string, #!# type
		:pdb_occupancy => :float,
		:pdb_temp_factor => :float,
		:pdb_segid => :string, #!# type
		:pdb_element => :string,
		:pdb_charge => :string, #!# type

		# MD
		:velocity => :float
	}

	# Properties: Hash class with check for allowed properties on write
	class Properties < Hash
		def []=(name, value)
			if Atom::ALLOWED_PROPERTIES.has_key?(name)
				super(name, value)
			else
				Cuby::error("Attempt to set undefined atom property ('#{name}')")
			end
		end
	end

	# Properties with special getter functions:
	# Atom.mass - if properties[:mass] is not set, returns default atomic mass from the periodic table

	#=======================================================================
	# Initialize
	#=======================================================================
	
	#> Atom.new(element, x = 0.0, y=0.0, z=0.0, properties = {}) # => Atom
	#!
	
	def initialize(element = :X, x = 0.0, y=0.0, z=0.0, properties = {})
		#: Atom constructor. 
		super(x,y,z)
		if element.class == Symbol
			Atom.check_element(element)
			@element = element
		elsif element.class == String
			@element = Atom.element_from_string(element)
		else
			raise(TypeError, "Element must be symbol or string")
		end
		@properties = Properties.new
		properties.each_pair {|p, value|
			@properties[p] = value
		}
	end
	#!!
	
	def Atom.from_coordinate(element, coordinate,  properties = {})
		return Atom.new(element, coordinate.x, coordinate.y, coordinate.z, properties)
	end
	
	def Atom.from_s(string)
		a = string.split
		Cuby::error("Error in string to atom conversion, less than 4 records found\n('#{string}')") if a.size < 4
		Cuby::error("Error in string to atom conversion, x coordinate is not a number\n('#{string}')") unless a[1] =~ /[-]?[0-9]+\.[0-9]+/
		Cuby::error("Error in string to atom conversion, y coordinate is not a number\n('#{string}')") unless a[2] =~ /[-]?[0-9]+\.[0-9]+/
		Cuby::error("Error in string to atom conversion, z coordinate is not a number\n('#{string}')") unless a[3] =~ /[-]?[0-9]+\.[0-9]+/
		at = Atom.new(a.shift, a.shift.to_f, a.shift.to_f, a.shift.to_f)
		# Properties
		a.each{|s|
			if ALLOWED_PROPERTIES.has_key?(s.downcase.to_sym)
				if ALLOWED_PROPERTIES[s.downcase.to_sym] == :boolean
					at.properties[s.downcase.to_sym] = true
				else
					Cuby::error("Error in string to atom conversion, property value not assigned '#{s}'") unless f.size == 2
				end
				next
			end

			if s =~ /=/
				f = s.split('=')
				Cuby::error("Error in string to atom conversion, incorrect property assignment '#{s}'") unless f.size == 2
				key =  f[0].downcase.to_sym
				if ALLOWED_PROPERTIES.has_key?(key)
					case ALLOWED_PROPERTIES[key]
					when :string
						at.properties[key] = f[1]
					when :float
						at.properties[key] = f[1].to_f
					when :integer
						at.properties[key] = f[1].to_i
					end

					# Convert mass to cuby units
					if  key == :mass
						at.properties[key] = at.properties[key] * GMOL2UNIT
					end

					next
				else
					Cuby::error("Error in string to atom conversion, unknown property '#{f[0]}'")
				end
			end

			Cuby::error("Error in string to atom conversion, unknown property '#{s}'")
		}
		return at
	end
	
	#!
	# disable Atom[] constructor
	def self.[](*array) # => Atom
		raise(NoMethodError,"This method is disabled")
	end
	#!!

	def self.from_array(element, array, properties = {}) # => Atom
		#: Constructor reading the coordinates from an array
		return Atom.new(element,array[0], array[1], array[2], properties)
	end

	def self.from_coordinate(element, coordinate, properties = {}) # => Atom
		#: Constructor reading the coordinates from Coordinate-like class
		return Atom.new(element,coordinate.x, coordinate.y, coordinate.z, properties)
	end

	def Atom.element_from_string(string) # => Symbol
		#: Method used to parse string valid element Symbol. Raises an exception
		#: when it is impossible to assign element from known periodic table.
		s = string.downcase.capitalize
		s = "X" if s == "Xx" # Two-letter dummy atom fix
		symbol = s.to_sym
		Atom.check_element(symbol)
		return symbol
	end

	def Atom.check_element(symbol) # => nil
		#: Checks whether string is valid element abbreviation, raises exception otherwise
		unless PeriodicTable::ELEMENTS.include?(symbol)
			if PeriodicTable::ELEMENTS.include?(symbol.to_s.downcase.capitalize.to_sym)
				raise("Wrong name of element - use capitalized form, i.e. :He for helium")
			else
				raise("Wrong name of element - unknown abbreviation #{symbol}")
			end
		end
	end

	#=======================================================================
	# Direct access to properties
	#=======================================================================
	
	def mass
		if @properties[:mass]
			return @properties[:mass]
		else
			return PeriodicTable::mass(@element)
		end
	end

	#=======================================================================
	# basics
	#=======================================================================
	
	def deep_copy # => Atom
		#: Real copy of the object and all its instance variables
		atom_copy = Atom.new(@element, self.x, self.y, self.z, @properties.deep_copy)
		return atom_copy
	end

	def dup # => Atom
		#: Real copy of the object and all its instance variables
		return deep_copy
	end

	def inspect
		return "Atom: #{to_s}"
	end

	#=======================================================================
	# type conversion
	#=======================================================================

	def to_s(format = "15.8f") # => String
		#: Converts the complete atom information to string
		props = ""
		@properties.each_pair{|name, value|
			# Original file index is not written
			next if name == :file_index

			# Mass is converted to g/mol
			if name == :mass
				value = value * UNIT2GMOL
			end

			if value.kind_of?(Numeric)
				props << " #{name}=#{value}"
			elsif value.class == String
				#!# escape characters if needed
				props << " #{name}=#{value}"
			elsif value.class == TrueClass
				props << " #{name}"
			else
				raise "Atom property printing in to_s method is limited to numbers, stringsi\nand boolean values,\nfound #{name} = #{value} (#{value.class})"
			end
		}
		return sprintf("%-4s%#{format}%#{format}%#{format}%s",@element.to_s, self.x, self.y, self.z, props)
	end

	#=======================================================================
	# comparison
	#=======================================================================
	
	def == (atom) # => true | false
		# Comparison operator, compares coordinates and element
		raise(TypeError,"Atom can be compared only to Atom, this is #{atom.class}") unless atom.kind_of?(Atom)
		return (self.x==atom.x && self.y==atom.y && self.z==atom.z && @element == atom.element)
	end

	def =~ (atom) # => true | false
		#: Loose comparison, difference between each pair of coordinate elements must be
		#: smaller than class variable epsilon and the element must be the same
		raise(TypeError,"Atom can be compared only to Atom") unless atom.kind_of?(Atom)
		return ((self.x-atom.x).abs < @@epsilon && 
		        (self.y-atom.y).abs < @@epsilon &&
		        (self.z-atom.z).abs < @@epsilon &&
		        @element == atom.element)
	end

	#!
	#=======================================================================
	# Methods removed
	#=======================================================================
	# Following inherited methods are disabled:

	def map!
		raise(NoMethodError,"This method is disabled")
	end
	#!!

	#=======================================================================
	# Element-related
	#=======================================================================
	
	def proton_number # => integer
		#: Returns proton number of the element. Dummy atom (symbol X) has
		#: proton number 0 here.
		return PeriodicTable::ELEMENTS.index(@element)
	end

	#=======================================================================
	# YAML support
	#=======================================================================
	
	require "yaml"
	begin
		yaml_engine = "psych" if YAML == Psych
		yaml_engine = "syck" if YAML == Syck
	rescue
		# ruby1.8 uses syck by default
		yaml_engine = 'syck'
	end

	if yaml_engine == 'syck'
		yaml_as "tag:cuby.molecular.cz,2009:#{self}"
		def to_yaml(opts = {})
			YAML::quick_emit( self.object_id, opts ) do |out|
				out.map(taguri) do |map|
					map.add('x', self.x)
					map.add('y', self.y)
					map.add('z', self.z)
					map.add('element', self.element)
					map.add('properties', self.properties) if self.properties.size > 0
				end 
			end
		end

		def self.yaml_new(klass, tag, val)
			if val['properties']
				Atom.new(val['element'], val['x'], val['y'], val['z'], val['properties'])
			else
				Atom.new(val['element'], val['x'], val['y'], val['z'])
			end
		end

	elsif yaml_engine == 'psych'
		# Psych yaml engine
	
		def encode_with(coder)
			coder.tag = "!cuby.molecular.cz,2009/Atom"
			coder.style = Psych::Nodes::Mapping::FLOW # BLOCK or FLOW
			coder['x'] = self.x
			coder['y'] = self.y
			coder['z'] = self.z
			coder['element'] = self.element
			if self.properties.size > 0
				coder['properties'] = self.properties
			end
		end

		def init_with(coder)
			# This is not used when custom tag is defined
		end

		YAML::add_domain_type('cuby.molecular.cz,2009', 'Atom') { |type, val|
			if val['properties']
				Atom.new(val['element'], val['x'], val['y'], val['z'], val['properties'])
			else
				Atom.new(val['element'], val['x'], val['y'], val['z'])
			end
		}
	else
		raise "YAML engine can not be determined"
	end

	#=======================================================================
	# Marshalling fix
	#=======================================================================
	# Marshalling seems to be broken when Atom inherits from binary version
	# of Coordinate (which itself serializes well). This code fixes that.
	
	def marshal_dump
		return [self.x, self.y, self.z, @element, @properties]
	end

	def marshal_load(array)
		x,y,z, @element, @properties = array
		self.x = x
		self.y = y
		self.z = z
	end
end
