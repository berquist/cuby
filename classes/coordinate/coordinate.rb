################################################################################
#
# class Coordinate - wrapper
#
# Author: Jan Rezac
# Date created: 2008-08-30
# License: Cuby license
# Description: Coordinate in 3D space
# Status: Documented, Tested
#
################################################################################

require "classes/cuby"

#: Basic (and thus efficient) implementation of 3D vector. This file loads either ruby version
#: or C extension, and adds methods from module CoordinateVector.

# Use binary if possible
if (FileTest.exists?(File.dirname(__FILE__) + "/coordinate_c.so") ||
   FileTest.exists?(File.dirname(__FILE__) + "/coordinate_c.bundle")) && !$coordinate_force_ruby
	# Ruby version match test
	compiled_version = IO.readlines(File.dirname(__FILE__) + "/_version")[0].strip
	if compiled_version != RUBY_VERSION
		raise "Extension 'coordinate' had been compiled with another version of ruby"
	end
	# Load the C version
	require "classes/coordinate/coordinate_c"
	Cuby.log.puts_debug("Coordinate class uses binary library")
else
	# Load ruby version
	require "classes/coordinate/coordinate_ruby"
	Cuby.log.puts_debug("Coordinate class uses ruby implementation")
end

class Coordinate
	#=======================================================================
	# Include CoordinateVector into Coordinate
	#=======================================================================
	require "classes/coordinate/coordinate_vector"
	include CoordinateVector

	#=======================================================================
	# Class methods - measurements
	#=======================================================================

	def Coordinate.angle(a,b,c)
		u = a - b
		v = c - b
		return u.angle(v)
	end

	def Coordinate.torsion(a,b,c,d)
		v1 = b - a
		v2 = c - b
		v3 = d - c
		d = Math::atan2((v1*v2.abs).dot(v2.cross_product(v3)), (v1.cross_product(v2)).dot(v2.cross_product(v3)))
		d += 2.0 * Math::PI if d <= -Math::PI
		d -= 2.0 * Math::PI if d > Math::PI
		return d
	end

	#=======================================================================
	# Polar coordinates
	#=======================================================================
	
	def to_polar
		#: Conversion to polar coordinates. 
		#: Yields array [r, theta, phi]
		r = self.abs
		theta = Math::acos(self.z/r)
		phi = Math::atan2(self.y,self.x)
		return [r, theta, phi]
	end

	def Coordinate.from_polar(r, theta, phi)
		#: Constructor from polar coordinates
		return self.new(
			r * Math::sin(theta) * Math::cos(phi),
			r * Math::sin(theta) * Math::sin(phi),
			r * Math::cos(theta)
		)
	end

	def Coordinate.from_polar_deg(r, theta, phi)
		#: Constructor from polar coordinates in degrees
		return Coordinate.from_polar(r, Math::PI * theta / 180, Math::PI * phi / 180)
	end

	#=======================================================================
	# Advanced measurements
	#=======================================================================
	
	def point_plane_distance(plane_normal_vec, plane_point)
		#: Calculate perprndicular distance of the point to a plane
		#: defined by its normal vector and one point that lies
		#: on the plane.
		return ((self - plane_point).dot(plane_normal_vec.normalize)).abs
	end

	def point_line_distance(line_point_1, line_point_2)
		#: Calculate perprndicular distance of the point to a line
		#: defined by two points
		return ((self - line_point_1).cross_product(self - line_point_2)).abs / 
			(line_point_2 - line_point_1).abs
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
				end
			end
		end

		def self.yaml_new(klass, tag, val)
			Coordinate.new(val['x'], val['y'], val['z'])
		end

	elsif yaml_engine == 'psych'
		# Psych yaml engine
	
		def encode_with(coder)
			coder.tag = "!cuby.molecular.cz,2009/Coordinate"
			coder.style = Psych::Nodes::Mapping::FLOW # BLOCK or FLOW
			coder['x'] = self.x
			coder['y'] = self.y
			coder['z'] = self.z
		end

		def init_with(coder)
			# This is not used when custom tag is defined
			self.x = coder['x']
			self.y = coder['y']
			self.z = coder['z']
		end

		YAML::add_domain_type('cuby.molecular.cz,2009', 'Coordinate') { |type, value|
			Coordinate.new(value['x'], value['y'], value['z'])
		}
	else
		raise "YAML engine can not be determined"
	end


end

################################################################################
#: Extension of numeric types to allow commutative multiplication
################################################################################

#===============================================================================
# Float class extension
#===============================================================================
class Float
	alias mul_without_coordinate *

	def *(operand)
		#: Multiplication is extended to handle Float * Coordinate
		if operand.kind_of?(Coordinate)
			return operand * self
		else
			return mul_without_coordinate(operand)
		end
	end
end

#===============================================================================
# Fixnum class extension
#===============================================================================
class Fixnum
	alias mul_without_coordinate *

	def *(operand)
		#: Multiplication is extended to handle Fixnum * Coordinate
		if operand.kind_of?(Coordinate)
			return operand * self
		else
			return self.mul_without_coordinate(operand)
		end
	end
end

