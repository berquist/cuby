# Series of geometries as read from the settings
# It expects 'settings' keyword to point to and .xyz file
# Template from 'geometry_template' is used when this keyword is set
# Labels are read from comment line of the .xyz file


class GeometrySeries < Array
	attr_reader :labels

	def initialize(settings, label_text="scan:")
		# Read template
		if settings.set?(:geometry_template)
			template = Geometry.from_file(settings[:geometry_template])
		end

		# Read geometries
		@labels = []
		f = File.open(File.expand_path(settings[:geometry]))
		g = Geometry.new
		while g.read_xyz(:file => f)
			if settings.set?(:geometry_template)
				g2 = template.deep_copy
				g2.copy_coordinates_from(g)
				self << g2
			else
				self << g
			end
			# Get label from second line of the xyz file
			label = g.info[:xyz_remark]
			if label =~ /^#{label_text}/
				@labels << label.gsub(/^#{label_text}\s*/,"")
			end
			g = Geometry.new
		end
		f.close

		if @labels.size > 0 # at least some labels found
			if @labels.size != size
				# When the number of labels does not correspond with a number of geometries,
				# remove them
				@labels = []
			end
		end
	end
end
