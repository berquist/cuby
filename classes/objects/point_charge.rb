################################################################################
#
# Class PointCharge
#
# Author: Jan Rezac
# Date created: 2008-09-05
# License: Cuby license
# Description: Point charge is a coordinate + charge
# Status: Works, Documented
#
################################################################################

#: Ancestor of Coordinate, adding attribute charge

require "classes/coordinate/coordinate.rb"

class PointCharge < Coordinate

	# List of privately used instance variables
	# @x, @y, @z - inherited from coordinate
	
	# List of public attributes
	attr_reader 	:charge

	# Additional information (optional, not initialized)

	# Backreference to the atom from which the charge was created
	attr_accessor	:from_atom

	#=======================================================================
	# Initialize
	#=======================================================================
	
	#> PointCharge.new(x = 0.0, y=0.0, z=0.0, charge = 0.0) # => PointCharge
	#!
	def initialize(x = 0.0, y=0.0, z=0.0, charge = 0.0)
		#: Constructor
		super(x,y,z)
		@charge = charge.to_f
	end
	#!!
	
	#!
	# disable PointCharge[] constructor
	def self.[](*array)
		raise(NoMethodError,"This method is disabled")
	end
	#!!

	def self.from_array(array,c = 0.0)
		#: Creates point charge from coordinates in an array and the
		#: charge.
		return PointCharge.new(array[0], array[1], array[2], c)
	end

	def self.from_coordinate(coordinate, c = 0.0)
		#: Creates point charge from coordinates in Coordinate object and the
		#: charge.
		return PointCharge.new(coordinate.x, coordinate.y, coordinate.z, c)
	end

	def self.from_s(string)
		#: Creates PointCharge from a string consisting of four numbers
		#: separated by whitespace, ordered as x,y,z,c

		# Check format
		a = string.strip.split
		Cuby::error("Point charge in format 'x y z charge' expected but less than 4 records found:\n#{string}") if a.size < 4
		Cuby::error("Point charge in format 'x y z charge' expected but more than 4 records found:\n#{string}") if a.size > 4
		a.each_index{|i|
			unless a[i] =~ /[-]?[0-9]+\.?[0-9]*/
				Cuby::error("Reading poin charge from string: record that doesn't look like a number found:\n#{a[i]}")
			end
		}

		# Return new PointCharge
		return PointCharge.new(a[0].to_f, a[1].to_f, a[2].to_f, a[3].to_f)
	end

	#=======================================================================
	# basics
	#=======================================================================
	
	def charge=(value)
		@charge = value.to_f
	end
	
	def deep_copy
		#: Deep copy of the object
		return PointCharge.new(self.x, self.y, self.z, @charge)
	end

	#=======================================================================
	# Operators
	#=======================================================================
	
	undef_method(:"-@")
	
	#=======================================================================
	# type conversion
	#=======================================================================

	def to_s(short = false)
		#: Conversion to String
		if short
			format = "%15.8f"
			return sprintf(format * 3 + "%10.4f", self.x, self.y, self.z, @charge)
		else
			return super() + ", charge #{@charge}"
		end

	end

	def inspect
		#: Use to_s to print nice output on inspect calls
		return to_s
	end

	def to_point_charge
		#: Allows conversion of descendants to PointCharge class
		return PointCharge.new(self.x,self.y,self.z,@charge)
	end

	#=======================================================================
	# comparison
	#=======================================================================
	
	def == (point_charge) # true | false
		#: Comparison operator, compares coordinates and charge
		raise(TypeError,"PointCharge can be compared only to PointCharge") unless point_charge.kind_of?(PointCharge)
		return (self.x==point_charge.x && self.y==point_charge.y && self.z==point_charge.z && @charge == point_charge.charge)
	end

	def =~ (point_charge) # true | false
		#: Loose comparison, difference between each pair of coordinate elements must be
		#: smaller than class variable epsilon and the charge must match exactly.
		raise(TypeError,"PointCharge can be compared only to PointCharge") unless point_charge.kind_of?(PointCharge)
		return ((self.x-point_charge.x).abs < @@epsilon && 
		        (self.y-point_charge.y).abs < @@epsilon &&
		        (self.z-point_charge.z).abs < @@epsilon &&
		        @charge == point_charge.charge)
	end

	#!
	#=======================================================================
	# Methods removed
	#=======================================================================
	
	def map!
		raise(NoMethodError,"This method is disabled")
	end
	#!!
	
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
					map.add('charge', self.charge)
					map.add('from_atom', self.from_atom) if self.from_atom
				end
			end
		end

		def self.yaml_new(klass, tag, val)
			c = PointCharge.new(val['x'], val['y'], val['z'], val['charge'])
			c.from_atom = val['from_atom'] if val['from_atom']
			c
		end

	elsif yaml_engine == 'psych'
		# Psych yaml engine
	
		def encode_with(coder)
			coder.tag = "!cuby.molecular.cz,2009/PointCharge"
			coder.style = Psych::Nodes::Mapping::FLOW # BLOCK or FLOW
			coder['x'] = self.x
			coder['y'] = self.y
			coder['z'] = self.z
			coder['charge'] = self.charge
			coder['from_atom'] = self.from_atom if self.from_atom
		end

		YAML::add_domain_type('cuby.molecular.cz,2009', 'PointCharge') { |type, val|
			c = PointCharge.new(val['x'], val['y'], val['z'], val['charge'])
			c.from_atom = val['from_atom'] if val['from_atom']
			c
		}
	else
		raise "YAML engine can not be determined"
	end
	
end
