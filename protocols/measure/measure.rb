require "classes/internal_coordinates/internal_coordinate"

module ProtocolMeasure
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :tools
	PROTOCOL_DESCRIPTION = "Performs measurements on the molecular geometry"
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Geometry measurements')
	end

	def prepare_protocol
		# Check whether geometry is an .xyz file
		unless FileTest.file?(File.expand_path(@settings[:geometry]))
			Cuby::error("For the measure job, geometry keyword must reference to a .xyz file")
		end
	end

	def run_protocol(queue)
	end

	def print_protocol
		# Get template geometry
		template = template_geometry(@settings)

		# Build measurements
		measurements = []
		@settings[:measurements].each_pair{|name, expression|
			m = MeasureExpr.from_string(expression, template)
			m.name = name
			measurements << m
		}

		# Print table header
		header = ''
		measurements.each {|e| 
			header += sprintf("%-12s",e.name) 
		}
		Cuby::log.puts_v(:normal, header)

		# Iterate over frames loaded into the template geometry
		iterate_frames(@settings, template) {|geometry|
			s = ''
			measurements.each {|e| s += sprintf("%-12.3f",e.calculate(geometry)) }
			Cuby::log.puts(s)
		}
	end

	def cleanup_protocol
	end

	def results
		return nil
	end


	#=======================================================================
	# Helpers
	#=======================================================================
	
	def template_geometry(settings)
		# Load the first geometry or a template
		geometry = Geometry.new
		if settings.set?(:geometry_template)
			geometry.read_file(File.expand_path(settings[:geometry_template]))
		else
			f = File.open(File.expand_path(settings[:geometry]))
			geometry.read_xyz(:file => f)
			f.close
		end
		return geometry
	end
	
	def iterate_frames(settings, template_geometry)
		# Iterate over the geometries
		f = File.open(File.expand_path(settings[:geometry]))
		g = Geometry.new
		while g.read_xyz(:file => f)
			template_geometry.copy_coordinates_from(g)
			yield(template_geometry)
		end
		return nil
	end
end

################################################################################
# The measurement object
################################################################################

class MeasureExpr < Struct.new(:name, :type, :atomlists, :unit, :calculator)
	def MeasureExpr.from_string(s, geometry)
		# look for long form with name
		match = /^\s*(\S+)\s+=\s+([a-z,A-Z,_]+)\s*\((.*)\)/.match(s)
		if match != nil
			name = match[1]
			type = match[2]
			params = match[3]
		else
			# try short form without name
			match = /^\s*([a-z,A-Z,_]+)\s*\((.*)\)/.match(s)
			if match != nil
				name = match[1]
				type = match[1]
				params = match[2]
			else
				Cuby::error "Wrong expression syntax in '#{s}'" if match == nil
			end
		end

		# Get type
		case type
		when /^dist/i
			typesymbol = :distance
			centers = 2
			unit = "A"
		when /^angle/i
			typesymbol = :angle
			centers = 3
			unit = "deg"
		when /^tor/i
			typesymbol = :torsion
			centers = 4
			unit = "deg"
		when /^closest_dist/i
			typesymbol = :closest_dist
			centers = 2
			unit = "A"
		when /^point_plane/i
			typesymbol = :point_plane
			centers = 4
			unit = "A"
		when /^point_normal/i
			typesymbol = :point_normal
			centers = 4
			unit = "A"
		when /^difference_d/i
			typesymbol = :distance_difference
			centers = 4
			unit = "A"
		else
			Cuby::error "Unknown measurement type '#{type}'"
		end

		# Parse atomlists
		selections = params.split(/\s*;\s*/) # Get selections
		Cuby::error "Wrong number of parameters for #{typesymbol.to_s}: #{selections.size} instead of #{centers}" unless centers == selections.size
		atomlists = selections.map{|sel| geometry.atomlist_from_selection(sel) }

		# Create calculator object if it is possible
		calculator = nil
		if [:distance, :angle, :torsion, :distance_difference].include?(typesymbol)
			calculator = GroupsInternalCoord.new(typesymbol, geometry, *atomlists)
		end

		# Create and return instance of MeasureExpr
		e = MeasureExpr.new(name, typesymbol, atomlists, unit, calculator)
		return e
	end

	def calculate(geometry)
		case self.type
		when :distance
			return self.calculator.value
		when :angle
			return self.calculator.value / Math::PI * 180.0
		when :torsion
			return self.calculator.value / Math::PI * 180.0
		when :closest_dist
			min, i, j = geometry.closest_contact(self.atomlists[0], self.atomlists[1])
			return min
		when :distance_difference
			return self.calculator.value
		when :point_plane
			p = geometry.center_of_mass_atomlist(self.atomlists[0])
			a = geometry.center_of_mass_atomlist(self.atomlists[1])
			b = geometry.center_of_mass_atomlist(self.atomlists[2])
			c = geometry.center_of_mass_atomlist(self.atomlists[3])
			vec = (a-b).cross_product(a-c)
			return p.point_plane_distance(vec, a)
		when :point_normal
			p = geometry.center_of_mass_atomlist(self.atomlists[0])
			a = geometry.center_of_mass_atomlist(self.atomlists[1])
			b = geometry.center_of_mass_atomlist(self.atomlists[2])
			c = geometry.center_of_mass_atomlist(self.atomlists[3])
			vec = (a-b).cross_product(a-c)
			center = a
			return p.point_line_distance(center, center + vec)
		end
	end
end
