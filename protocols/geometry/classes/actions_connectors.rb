module ActionsConnectors

	#=======================================================================
	# Actions - adding connectors
	#=======================================================================

	def add_connector
		# Add a connector dummy atom. It's position is defined by connector_orientation keyword
		case @settings[:connector_type]
		when :linear
			sel_no = 3
		when :perpendicular
			sel_no = 3
		when :apex
			sel_no = 1
		end

		# Separate selections
		selections = @settings[:connector_orientation].strip.split(";")
		unless selections.size == sel_no
			Cuby::error "Keyword connector_orientation must define #{sel_no} selections"
		end

		# Convert them to atomlists
		atomlists = selections.map{|s|
			list = @geometry.atomlist_from_selection(s)
			Cuby::error "Selection in connector_orientation selects no atoms" unless list.size > 0
			list
		}


		# Minimum distance 0.6
		conn_dist = 0.6
		# Or largest vdw radius in the selection
		atomlists[0].each{|i|
			r = PeriodicTable.vdw_radius(@geometry[i].element)
			if r > conn_dist
				conn_dist = r
			end
		}
		case @settings[:connector_type]
		when :linear
			@geometry.orient!(atomlists[0],atomlists[1],atomlists[2])
			@geometry << Atom.new(:X, conn_dist - 0.5, 0, 0)
			@geometry << Atom.new(:X, conn_dist, 0, 0)
		when :perpendicular
			@geometry.orient!(atomlists[0],atomlists[1],atomlists[2])
			@geometry << Atom.new(:X, 0, 0, conn_dist - 0.5)
			@geometry << Atom.new(:X, 0, 0, conn_dist)
		when :apex
			unless atomlists[0].size == 1
				Cuby::error "Connctor of type 'apex' is defined only for a single atom"
			end
			atom_i = atomlists[0].first
			conn = Connectivity.new(@geometry)
			vec = Coordinate[0,0,0]
			conn.bound_atoms(atom_i).each{|atom_j|
				bond = @geometry[atom_i] - @geometry[atom_j]
				vec += bond / bond.abs
			}
			# Normalize the vector
			vec = vec / vec.abs
			@geometry << Atom.from_coordinate(:X, @geometry[atom_i] + vec * (conn_dist-0.5))
			@geometry << Atom.from_coordinate(:X, @geometry[atom_i] + vec * conn_dist)
		end

		# The connector is the last atom - write it to header
		@geometry.info[:protocol_comment] = "connector=#{@geometry.size}"

	end

	#=======================================================================
	# Actions - structure building using connectors
	#=======================================================================
	
	def connect_connectors
		geometry2 = Geometry.from_file(@settings[:geometry2])

		# Get connectors
		if @settings[:geometry_connector].downcase == "auto"
			if @geometry.info[:xyz_remark] && (matchdata = /connector *= *(\S+)/.match(@geometry.info[:xyz_remark]))
				puts matchdata[1]
				@settings[:geometry_connector] = matchdata[1]
			else
				Cuby::error "Connector information not found in geometry"
			end
		end
		conn1 = @geometry.atomlist_from_selection(@settings[:geometry_connector])
		if @settings[:geometry_connector2].downcase == "auto"
			if geometry2.info[:xyz_remark] && (matchdata = /connector *= *(\S+)/.match(geometry2.info[:xyz_remark]))
				puts matchdata[1]
				@settings[:geometry_connector2] = matchdata[1]
			else
				Cuby::error "Connector information not found in geometry2"
			end
		end
		conn2 = geometry2.atomlist_from_selection(@settings[:geometry_connector2])
		Cuby::error("Keyword geometry_connector should select one atom (#{conn1.size} atoms selected)") unless conn1.size == 1
		Cuby::error("Keyword geometry_connector2 should select one atom (#{conn2.size} atoms selected)") unless conn2.size == 1
		c1_atom = @geometry[conn1[0]]
		c2_atom = geometry2[conn2[0]]

		puts "Connectors:"
		puts c1_atom
		puts c2_atom

		# Get anchors - closest atom to connector
		a1_atom = connect_connectors_get_anchor(@geometry, c1_atom)
		a2_atom = connect_connectors_get_anchor(geometry2, c2_atom)

		puts "Anchors:"
		puts a1_atom
		puts a2_atom
		puts

		# Move molecule 1 connector to origin
		@geometry.translate!(-(c1_atom.to_coordinate))

		# Move molecule 2
		geometry2.translate!(-(c2_atom.to_coordinate))

		# Align a1 c1=c2 a2 to line
		while Coordinate.angle(a1_atom, c1_atom, a2_atom) < 0.1
			# break linearity
			geometry2.rotate_angle_axis!(0.3, Coordinate[rand,rand,rand])
		end
		normal = (a1_atom - c1_atom).cross_product(a2_atom - c2_atom)
		normal = normal / normal.abs
		angle = Coordinate.angle(a1_atom, c1_atom, a2_atom)
		geometry2.rotate_angle_axis!((Math::PI - angle), normal)
		angle = Coordinate.angle(a1_atom, c1_atom, a2_atom)

		# Rotate molecules around the axis
		com1 = @geometry.center
		com2 = geometry2.center
		torsion = Coordinate.torsion(com1,a1_atom,a2_atom,com2)
		# Maximize overlap - COMs must be aligned
		geometry2.rotate_angle_axis!(torsion, (a1_atom - a2_atom))

		if @settings[:geometry_connect_orientation] == :max_center_dist
			# Minimize overlap - COMs at opposite sides
			geometry2.rotate_angle_axis!(Math::PI, (a1_atom - a2_atom))
		end


		# Perform distance scan if requested
		if @settings[:geometry_connect_scan_vdw].size > 0
			scan_factors = @settings[:geometry_connect_scan_vdw]

			# Intramolecular axis - positive direction = increase distance
			intermolecular_v = (a2_atom - a1_atom)
			intermolecular_v /= intermolecular_v.abs

			# Geometries without dummy atoms
			g1 = Geometry.new
			@geometry.each{|atom|
				g1 << atom unless atom.element == :X
			}
			g2 = Geometry.new # Use deep copy of atoms here
			geometry2.each{|atom|
				g2 << atom.deep_copy unless atom.element == :X
			}
			# Move m2 to vdw contact
			closest_contact = 9.9e99
			ai = nil; aj = nil
			g1.each_index{|i|
				g2.each_index{|j|
					if g1[i].distance(g2[j]) < closest_contact
						closest_contact = g1[i].distance(g2[j])
						ai = i
						aj = j
					end
				}
			}
			vdw_distance =  PeriodicTable.vdw_radius(g1[ai].element) + PeriodicTable.vdw_radius(g2[aj].element)
			g2.translate!(intermolecular_v * (vdw_distance - closest_contact))

			g2.each{|atom|
				g1 << atom
			}
			scan_factors.each{|scale|
				translate = -(1.0 - scale) * vdw_distance
				g2.translate!(intermolecular_v * translate)
				g1.write_xyz(:file => "scan.xyz", :append => (scale != scan_factors.first))
				#g1.write_xyz(:file => "point_#{'%.2f' % scale}.xyz")
				g2.translate!(intermolecular_v * -translate)
			}
		end

		# Add geometry2 to geometry
		geometry2.each{|atom|
			@geometry << atom
		}
	end

	# Not an action, used in the method above
	def connect_connectors_get_anchor(geometry, atom_i) 
		a = nil
		mindist = 9.9e99
		geometry.each_index{|i|
			next if geometry[i] == atom_i
			next if geometry[i].element == :X
			if (d = atom_i.distance(geometry[i])) < mindist
				mindist = d
				a = i
			end
		}

		# If there are multiple atoms close, make an average
		atoms = []
		allowance = 1.1
		geometry.each_index{|i|
			next if geometry[i] == atom_i
			next if geometry[i].element == :X
			d = atom_i.distance(geometry[i])
			if d < mindist * allowance
				atoms << geometry[i]
			end
		}

		# Only one atom found - return it
		return atoms[0] if atoms.size == 1

		# Make average
		crd = Coordinate[0,0,0]
		atoms.each{|x|
			crd += x / atoms.size
		}
		a = Atom.from_coordinate(:X, crd)
		geometry << a
		return a
	end
end
