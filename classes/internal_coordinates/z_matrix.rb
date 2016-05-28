class ZMatrix
	class ZMatrixError < RuntimeError
	end

	class Coord

		attr_reader :ref_atom
		attr_accessor :value
		attr_accessor :name
		attr_reader :type


		def initialize(ref_atom, value_or_name, type)
			@type = type

			if value_or_name.kind_of?(Numeric)
				@value = value_or_name.to_f
				@value = @value / 180.0 * Math::PI unless @type == :distance

				@name = nil
			elsif value_or_name.kind_of?(String)
				@value = nil
				@name = value_or_name
			else
				raise(ZMatrixError, "Argument value_or_name must be either string or a number")
			end

			@ref_atom = ref_atom
		end

		def to_s
			s = sprintf("  %-4d", @ref_atom + 1)
			if @name
				s += sprintf("%8s", @name)
			else
				# Print z-matrix with angles in degrees
				val = @value
				val = val * 180.0 / Math::PI unless @type == :distance
				s += sprintf("%8.3f", val)
			end
			return s
		end

	end

	class Variable
		attr_accessor :value
		attr_accessor :type

		def initialize(value, type)
			@value, @type = value, type
		end
	end

	class Line
		attr_reader :element, :coordinates
		attr_reader :atom_in_geo

		def initialize(element, coordinates = [], atom_in_geo = nil)
			@element = element
			@coordinates = coordinates
			@atom_in_geo = atom_in_geo
		end

		def Line.build_deg(element, *args)
			case args.size
			when 0
				return Line.new(element)
			when 2
				return Line.new(element, [Coord.new(args[0], args[1], :distance)])
			when 4
				return Line.new(element, [Coord.new(args[0], args[1], :distance), Coord.new(args[2], args[3], :angle)])
			when 6
				return Line.new(element, [Coord.new(args[0], args[1], :distance), Coord.new(args[2], args[3], :angle), Coord.new(args[4], args[5], :torsion)])
			else
				raise(ZMatrixError, "Use 0, 2, 4 or 6 arguments to construct a z-matrix")
			end
		end

		def Line.from_s(string)
			words = string.strip.split
			element = Atom.element_from_string(words.shift)
			unless [0,2,4,6].include?(words.size)
				raise(ZMatrixError, "Wrong number of records in a z-matrix line '#{string}'")
			end
			args = []
			(words.size / 2).times{|i|
				refatom = words[i*2]
				value = words[i*2 + 1]

				# The reference atom
				if refatom =~ /^[0-9]+/
					refatom = refatom.to_i
					if refatom == 0
						Cuby::error("Numbering of atoms in z-matrix starts from 1 but zero was found")
					end
					refatom = refatom - 1
				else
					Cuby::error("Reference to a previous atom in z-matrix must be a number")
				end

				# Value or variable name
				if value =~ /^[-]?[0-9]*\.?[0-9]+$/
					value = value.to_f
				end

				args << refatom
				args << value
			}
			return Line.build_deg(element, *args)
		end

		def size
			@coordinates.size
		end

		def to_s
			s = sprintf("%-4s", @element.to_s)
			@coordinates.each{|crd| s += crd.to_s}
			return s
		end

		def distance
			return @coordinates[0]
		end

		def angle
			return @coordinates[1]
		end

		def torsion
			return @coordinates[2]
		end

	end

	#=======================================================================
	# ZMatrix class itself
	#=======================================================================

	attr_reader :lines

	def initialize
		@lines = []
	end

	def size
		return @lines.size
	end

	def add_line(line)
		case @lines.size
		when 0
			raise(ZMatrixError, "First line of z-matrix should contain no coordinates") unless line.size == 0
		when 1
			raise(ZMatrixError, "Second line of z-matrix should contain only distance") unless line.size == 1
		when 2
			raise(ZMatrixError, "Second line of z-matrix should contain only distance and angle") unless line.size == 2
		else
			raise(ZMatrixError, "A line of z-matrix should contain distance, angle and torsion coordinates") unless line.size == 3
		end
		@lines << line
	end

	def each_coord
		@lines.each{|line|
			line.coordinates.each{|coord|
				yield coord
			}
		}
	end

	def each_coord_with_index
		i = 0
		@lines.each{|line|
			line.coordinates.each{|coord|
				yield coord, i
				i += 1
			}
		}
	end

	def coord_by_index(i)
		case i
		when 0
			return @lines[1].coordinates[0]
		when 1
			return @lines[2].coordinates[0]
		when 2
			return @lines[2].coordinates[1]
		else
			fi = i - 3
			line = fi / 3 + 3
			return @lines[line].coordinates[fi % 3]
		end
	end

	def bond_list
		# Bond list as array of [i,j] arrays
		bondlist = []
		(@lines.size - 1).times {|i|
			bondlist << [i+1, @lines[i+1].distance.ref_atom]
		}
		return bondlist
	end

	def list_of_variables
		variables = {}
		each_coord{|coord|
			if coord.name
				if variables.has_key?(coord.name)
					unless coord.value == variables[coord.name].value
						raise "Two z-matrix coordinates with the same name have different values"
					end
				else
					variables[coord.name] = Variable.new(coord.value, coord.type)
				end
			end
		}
		return variables
	end

	def build_variables!
		count_d = 1
		count_a = 1
		count_t = 1
		each_coord{|coord|
			next if coord.name
			case coord.type
			when :distance
				coord.name = "d#{count_d}"
				count_d += 1
			when :angle
				coord.name = "a#{count_a}"
				count_a += 1
			when :torsion
				coord.name = "t#{count_t}"
				count_t += 1
			end
		}
		return nil
	end

	def to_s
		s = []
		@lines.each{|line|
			s << line.to_s
		}
		s = s.join("\n")

		# Collect names of variables
		variables = list_of_variables

		# Print variable names if there are any
		if variables.size > 0
			s += "\n\n"
			variables.each_pair{|name, var|
				value = var.value
				value = value * 180.0 / Math::PI unless var.type == :distance
				s += sprintf("%-12s%12.5f\n",name,value)
			}
		end
		return s
	end

	def ZMatrix.from_s(s)
		zmat = ZMatrix.new
		variables = {} # here, Variable class does not have to be used, we just collect the values

		mode = :atoms
		s.split("\n").each{|line|
			if line =~ /^\s*$/
				mode = :vars unless zmat.size == 0
				next
			end

			case mode
			when :atoms
				zmat.add_line(Line.from_s(line))
			when :vars
				name, value = line.split
				name.strip!
				value = value.to_f
				variables[name] = value
			end
		}

		# Assign values from variables list to coordinates
		zmat.each_coord{|coord|
			if coord.name
				unless variables[coord.name]
					Cuby::error "Z-matrix variable '#{coord.name}' does not have a value defined"
				end
				coord.value = variables[coord.name]
				# Angles are in degrees, convert them:
				coord.value = coord.value / 180.0 * Math::PI unless coord.type == :distance
			end
		}

		return zmat

	end

	def to_geometry
		geom = Geometry.new
		@lines.each_index{|i|
			line = @lines[i]
			case i
			when 0
				# just atom
				geom << Atom.new(line.element,0.0,0.0,0.0)
			when 1
				# atom in a distance
				geom << Atom.new(line.element,line.distance.value,0.0,0.0)
			when 2
				# distance & angle
				coord1 = geom[line.distance.ref_atom].to_coordinate
				coord2 = geom[line.angle.ref_atom].to_coordinate
				vec1 = coord1 - coord2 # bond vector
				axis = Coordinate.new(0.0, 0.0, 1.0) # rotate vector around z axis
				vec2 = vec1.rotate_angle_axis(Math::PI - line.angle.value,axis).normalize
				newcoord = vec2 * line.distance.value + coord1
				geom << Atom.new(line.element, newcoord.x, newcoord.y, newcoord.z)
			else
				# full specification
				coord1 = geom[line.distance.ref_atom].to_coordinate
				coord2 = geom[line.angle.ref_atom].to_coordinate
				coord3 = geom[line.torsion.ref_atom].to_coordinate
				vec1 = coord1 - coord2 # bond vector
				vec2 = coord2 - coord3 # torsion bond vector
				orig_plane_norm = vec1.cross_product(vec2) # plane of the three reference atoms
				axis = orig_plane_norm.rotate_angle_axis(line.torsion.value - Math::PI,vec1).normalize # apply tosrion angle to original plane norm
				vecnew = vec1.rotate_angle_axis(Math::PI - line.angle.value,axis).normalize # apply angle
				newcoord = vecnew * line.distance.value + coord1
				geom << Atom.new(line.element, newcoord.x, newcoord.y, newcoord.z)
			end
		}
		return geom
	end

	def update_geometry!(geometry)
		# Check size
		raise "Update of geometry from z-matrix: numbers of atoms are different" unless geometry.size == @lines.size

		@lines.each_index{|i|
			line = @lines[i]

			# Check element
			raise "Element mismatch on update of geometry from z-matrix" unless geometry[i].element == line.element

			case i
			when 0
				# just atom
				geometry[i].set_coord(Coordinate[0.0,0.0,0.0])
			when 1
				# atom in a distance
				geometry[i].set_coord(Coordinate[line.distance.value,0.0,0.0])
			when 2
				# distance & angle
				coord1 = geometry[line.distance.ref_atom].to_coordinate
				coord2 = geometry[line.angle.ref_atom].to_coordinate
				vec1 = coord1 - coord2 # bond vector
				axis = Coordinate.new(0.0, 0.0, 1.0) # rotate vector around z axis
				vec2 = vec1.rotate_angle_axis(Math::PI - line.angle.value,axis).normalize
				newcoord = vec2 * line.distance.value + coord1
				geometry[i].set_coord(newcoord)
			else
				# full specification
				coord1 = geometry[line.distance.ref_atom].to_coordinate
				coord2 = geometry[line.angle.ref_atom].to_coordinate
				coord3 = geometry[line.torsion.ref_atom].to_coordinate
				vec1 = coord1 - coord2 # bond vector
				vec2 = coord2 - coord3 # torsion bond vector
				orig_plane_norm = vec1.cross_product(vec2) # plane of the three reference atoms
				axis = orig_plane_norm.rotate_angle_axis(line.torsion.value - Math::PI,vec1).normalize # apply tosrion angle to original plane norm
				vecnew = vec1.rotate_angle_axis(Math::PI - line.angle.value,axis).normalize # apply angle
				newcoord = vecnew * line.distance.value + coord1
				geometry[i].set_coord(newcoord)
			end
		}
		return nil
	end

	#=======================================================================
	# Mopac format
	#=======================================================================

	def to_s_mopac(label = :none)
		str = []
		@lines.each_index{|i|
			line = @lines[i]
			# C         1.45099200 +1    0.0000000 +0    0.0000000 +0     1     0     0
			case label
			when :number
				s = sprintf("  %-10s", line.element.to_s + "(#{line.atom_in_geo+1})")
			when :user_bond
				bondlabel = ""
				if line.distance && line.distance.name =~ /^bond_/
					bondlabel = "(" + line.distance.name + ")"
				end
				s = sprintf("  %-10s", line.element.to_s + bondlabel)
			else
				s = sprintf("  %-3s", line.element.to_s)
			end
			dat = []
			line.coordinates.each{|crd|
				# Numbering of atoms starts at 1
				value = crd.value
				value = value * 180.0 / Math::PI unless crd.type == :distance
				dat << [value, 1, crd.ref_atom + 1]
			}
			while dat.size < 3
				dat << [0.0, 0, 0]
			end
			# Add coordinate columns
			dat.each{|col|
				s += sprintf("%15.8f +%1d", col[0], col[1])
			}
			# At refernce atom columns
			dat.each{|col|
				s += sprintf("%6d", col[2])
			}
			str << s
		}
		str = str.join("\n")
		return str
	end

	#=======================================================================
	# Z-matrix building
	#=======================================================================

	LINEAR_ANGLE_THRESHOLD = 179.0

	def ZMatrix.from_geometry_keep_order(geometry)
		zmat = ZMatrix.new
		variables = {}

		connectivity = Connectivity.new(geometry)

		geometry.each_index{|i|
			atom = geometry[i]
			case i
			when 0
				zmat.add_line(Line.new(atom.element))
			when 1
				dist = Coord.new(i-1, atom.distance(geometry[i-1]), :distance)
				zmat.add_line(Line.new(atom.element, [dist]))
			when 2
				dist = Coord.new(i-1, atom.distance(geometry[i-1]), :distance)
				angle = Coord.new(i-2, Coordinate.angle(atom, geometry[i-1], geometry[i-2]) / Math::PI * 180, :angle)
				if angle.value > 179.0
					Cuby::error("Z-matrix generation: linear angle encountered")
				end
				zmat.add_line(Line.new(atom.element, [dist, angle]))
			else
				dist = Coord.new(i-1, atom.distance(geometry[i-1]), :distance)
				angle = Coord.new(i-2, Coordinate.angle(atom, geometry[i-1], geometry[i-2]) / Math::PI * 180, :angle)
				torsion = Coord.new(i-3, Coordinate.torsion(atom, geometry[i-1], geometry[i-2], geometry[i-3]) / Math::PI * 180, :torsion)
				zmat.add_line(Line.new(atom.element, [dist, angle, torsion]))
			end
		}

		return zmat

	end

	def ZMatrix.from_geometry(geometry, bondlist = [])
		zmat = ZMatrix.new

		# Build connectivity
		connectivity = Connectivity.new(geometry)

		# Add all bonds from the bond list to connectivity
		if bondlist.size > 0
			bondlist.each{|bond|
				i, j = bond
				connectivity.add_bond(i,j)
			}
		end

		# Connect all molecules by their closest contacts
		molecules = connectivity.molecules
		while molecules.size > 1
			dist, mol, ai, aj = geometry.find_closest_molecule(molecules,0)
			connectivity.add_bond(ai,aj)
			molecules = connectivity.molecules
		end

		# Find main and side chains
		if bondlist.size == 0
			# Find longest chain and rank the atoms by distance from one end of the chain
			graph_tree_ranks = ZMatrix.find_longest_chain(geometry, connectivity)
			# Build all paths, starting from the longest
			all_paths = ZMatrix.paths_in_tree(geometry, connectivity, graph_tree_ranks)
		else
			# Find a chain that goes through all specified bonds

			# Build a skeleton from all specified bonds and their shortest connections
			# Bond clusters
			bond_clusters = []
			bond_clusters_bondlists = []
			Cuby::log.puts_debug "Processing bond lists:"
			bondlist.each{|bond|
				# Is any of the atoms in a cluster already?
				connects_to = []
				bond_clusters.each_index{|ci|
					overlap = bond_clusters[ci]  & bond
					if overlap.size == 1
						connects_to << ci
					elsif overlap.size > 1
						# If it connect two atoms from the same cluster, it closes a ring
						Cuby::error "Selected bonds make a ring, z-matrix could not be built"
					end
				}
				Cuby::log.puts_debug "   #{bond[0] + 1},#{bond[1] + 1} connects to #{connects_to.size}"
				case connects_to.size
				when 0
					# Not connected to any cluster, make a new one
					bond_clusters << bond
					bond_clusters_bondlists << [bond]
				when 1
					# Add the bond to a cluster
					bond_clusters[connects_to[0]] = bond_clusters[connects_to[0]] | bond
					bond_clusters_bondlists[connects_to[0]] << bond
				when 2
					# Connects two clusters - merge them
					c1 = bond_clusters[connects_to[0]]
					c2 = bond_clusters[connects_to[1]]
					bond_clusters.delete(c1)
					bond_clusters.delete(c2)
					bond_clusters << (c1 | c2)
					# merge bond lists
					l1 = bond_clusters_bondlists[connects_to[0]]
					l2 = bond_clusters_bondlists[connects_to[1]]
					bond_clusters_bondlists.delete(l1)
					bond_clusters_bondlists.delete(l2)
					bond_clusters_bondlists << (l1 | l2 | [bond])
				else
					raise "This should not happen"
				end
			}

			# Debug: save cluster geometries
			#bond_clusters.each_index{|i|
			#	geometry.write_file("cluster_#{i}.mol2", :mol2, {:bondlist => bond_clusters_bondlists[i]})
			#}

			# Connect the clusters to make a skeleton
			# Make a graph
			graph = MatrixGraph.empty(geometry.size)
			geometry.size.times{|i|
				i.times{|j|
					graph.add_edge_undirected(i,j, 1) if connectivity.is_bond?(i,j)
				}
			}

			# Connect the clusters by shortest paths
			bondlist_flat = bondlist.flatten
			while bond_clusters.size > 1
				shortest_dist = 9999999999999
				shortest_path = nil
				connects = nil
				# Find shortest connection
				bond_clusters.each_index{|ci|
					ci.times{|cj|
						bond_clusters[ci].each{|at_i|
							bond_clusters[cj].each{|at_j|
								dist, path = graph.dijkstra_best_path(at_i, at_j)
								# Skip if path goes throug cluster
								next if (path & bond_clusters[ci]).size > 1
								next if (path & bond_clusters[cj]).size > 1
								# Save if it is shortest
								if dist < shortest_dist
									shortest_dist = dist
									shortest_path = path
									connects = [ci, cj]
								end
							}
						}
					}
				}
				unless connects
					raise "Can not connect two clusters of bonds"
				end
				# Connect the clusters
				c1 = bond_clusters[connects[0]]
				c2 = bond_clusters[connects[1]]
				bond_clusters.delete(c1)
				bond_clusters.delete(c2)
				bond_clusters << (c1 | c2 | shortest_path)
				# Add the connecting path
				l1 = bond_clusters_bondlists[connects[0]]
				l2 = bond_clusters_bondlists[connects[1]]
				bond_clusters_bondlists.delete(l1)
				bond_clusters_bondlists.delete(l2)
				path_bl = []
				(shortest_path.size-1).times{|i|
					path_bl << [shortest_path[i], shortest_path[i+1]]
				}
				bond_clusters_bondlists << (l1 | l2 | path_bl)
			end
			skeleton = bond_clusters.first
			skeleton_bonds = bond_clusters_bondlists.first
			geometry.write_file("skeleton.mol2", :mol2, {:bondlist => bond_clusters_bondlists.first})

			# Connectivity of the skeleton only
			connectivity_skeleton = Connectivity.new(geometry.size)
			skeleton_bonds.each{|bond|
				a,b = bond
				connectivity_skeleton.add_bond(a,b)
			}

			# Rank atoms in the skeleton - there are no rings so all bonds will appear in the z-matrix
			start = ZMatrix.outermost_terminal_atom(geometry, connectivity_skeleton)
			graph_tree_ranks = Array.new(geometry.size, nil)
			ZMatrix.flood_ranks(graph_tree_ranks, geometry, connectivity_skeleton, start, skeleton)
			# Build path(s) for the skeleton only
			skeleton_paths = ZMatrix.paths_in_tree(geometry, connectivity_skeleton, graph_tree_ranks)

			# Grow sidechains on the skeleton
			ZMatrix.flood_ranks_add(graph_tree_ranks, geometry, connectivity)
			# And add them to the paths
			all_paths = ZMatrix.paths_in_tree(geometry, connectivity, graph_tree_ranks, skeleton_paths)
		end


		# all_paths contain all paths from a root of a tree to the branches
		# the root has rank 0, branch with highest rank is longest
		# puts all_paths.map{|x| x.map{|y| y+1}.join(",")}

		# Build z-matrix
		available_atoms = (0..(geometry.size-1)).to_a
		all_paths.each{|path|
			path.each_index{|i|
				atom_index = path[i]
				# Add only available atoms
				next unless available_atoms.include?(atom_index)
				available_atoms.delete(atom_index)

				# Add current atom
				atom = geometry[atom_index]
				coords = []

				# Distance coordinate
				if zmat.size > 0
					# Easy, always use previous atom in the path
					raise "No distance reference" if i < 1
					refatom_dist = path[i-1]
					refatom_dist_zmi = zmat.find_atom_index(refatom_dist)
					coord = Coord.new(refatom_dist_zmi, atom.distance(geometry[refatom_dist]), :distance)
					if bondlist.include?([atom_index, refatom_dist]) || bondlist.include?([refatom_dist, atom_index])
						bi = bondlist.index([atom_index, refatom_dist])
						bi = bondlist.index([refatom_dist, atom_index]) unless bi
						coord.name = "bond_#{bi+1}"
					end
					coords << coord
				raise "Distance refatom not found" unless refatom_dist_zmi
				end

				# Angle coordinate
				if zmat.size > 1
					refatom_angle_zmi = nil
					refatom_angle = nil
					if i < 2
						# No angle reference, search for one
						neighbors = connectivity.bound_atoms(refatom_dist)
						neighbors.delete(atom_index)
						raise "No angle reference" if neighbors.size == 0
						neighbors.each{|n|
							refatom_angle = n
							refatom_angle_zmi = zmat.find_atom_index(refatom_angle)
							break if refatom_angle_zmi
						}
						raise "Angle refatom not found in zmat (outside path)" unless refatom_angle_zmi
					else
						refatom_angle = path[i-2]
						refatom_angle_zmi = zmat.find_atom_index(refatom_angle)
						raise "Angle refatom not in zmat (in path)" unless refatom_angle_zmi
					end

					# Check for near-linear angle
					angle = Coordinate.angle(atom, geometry[refatom_dist], geometry[refatom_angle]) / Math::PI * 180
					if angle > LINEAR_ANGLE_THRESHOLD
						# Close to linear, dumy atom have to be added
						# Bond vector
						vec = geometry[refatom_dist] - geometry[refatom_angle]
						vec = vec / vec.abs
						# Choose second vector not parallel with the bond
						if i < 3 # Only two previous atoms available
							vec2 = Coordinate[1,0,0]
							vec2 = Coordinate[0,1,0] if vec.angle(vec2) < 0.1
						else
							vec2 = geometry[refatom_angle] - geometry[path[i-3]]
							vec2 = Coordinate[1,0,0] if vec.angle(vec2) < 0.1
							vec2 = Coordinate[0,1,0] if vec.angle(vec2) < 0.1
						end
						# Find normal vector to the plane
						normal = vec.cross_product(vec2)
						normal = normal / normal.abs
						# Add atom perpendicular to the bond at the end of geometry
						atom_x = Atom.from_coordinate(:X, geometry[refatom_dist] + normal)
						atom_x_i = geometry.size
						geometry[atom_x_i] = atom_x
						# Add the dumy atom and bond to it bond to connectivity
						connectivity = Connectivity.resized(connectivity,1)
						connectivity.add_bond(refatom_dist, atom_x_i)

						# Add the dummy atom to the z-matrix
						coords_x = [
							Coord.new(zmat.find_atom_index(refatom_dist), 1.0, :distance),
							Coord.new(zmat.find_atom_index(refatom_angle), 90.0, :angle)

						]
						# Add a torsion for the dummy atom if needed
						if zmat.size > 2
							x_torsion_ref = nil
							# If there is atom available in the path, try using it
							if i >= 3
								x_torsion_ref = path[i-3] 
								if Coordinate.angle(geometry[refatom_dist], geometry[refatom_angle], geometry[x_torsion_ref]) / Math::PI * 180 > LINEAR_ANGLE_THRESHOLD
									x_torsion_ref = nil
									linear = true
								else
									unless zmat.find_atom_index(x_torsion_ref)
										raise "Ref not in zmat"
									end
								end
							end
							# If it is not there, or if it is in a linear angle, find another one
							unless x_torsion_ref
								neighbors = connectivity.bound_atoms(refatom_angle)
								neighbors.delete(refatom_dist)
								neighbors.delete(x_torsion_ref) if x_torsion_ref
								x_torsion_ref = nil
								neighbors.each{|ni|
									if zmat.find_atom_index(ni)
										x_torsion_ref = ni
										break
									end
								}
								unless x_torsion_ref
									raise "Building a torsion for dummy atom added to linear angle failed (zmat.size=#{zmat.size}, i=#{i}, linear=#{linear})"
								end

							end
							# Add the torsion to coordinates defining the dummy atom
							x_t_a = Coordinate.torsion(atom_x, geometry[refatom_dist], geometry[refatom_angle],  geometry[x_torsion_ref]) / Math::PI * 180
							refatom_t = zmat.find_atom_index(x_torsion_ref)
							raise "Dummy atom torsion not defined" unless refatom_t
							coords_x << Coord.new(refatom_t, x_t_a, :torsion)
						end

						# Add the dummy atom to z-matrix
						zmat.add_line(Line.new(atom_x.element, coords_x, atom_x_i))

						# Use the dummy atom to define the angle that was linear
						refatom_angle = atom_x_i
						refatom_angle_zmi = zmat.find_atom_index(atom_x_i)
						angle = angle - 90.0
					end

					# Add the angle
					coords << Coord.new(refatom_angle_zmi, angle, :angle)
				end

				# Torsion coordinate
				if zmat.size > 2
					refatom_torsion = nil
					# If there is atom available in the path, use it
					if i >= 3
						refatom_torsion = path[i-3]
						# If the atom is in a linear angle, discard it
						if Coordinate.angle(geometry[refatom_dist], geometry[refatom_angle], geometry[refatom_torsion]) / Math::PI * 180 > LINEAR_ANGLE_THRESHOLD
							refatom_torsion = nil
							was_linear = true
						end
					end

					# Look for another atom bound to the one defining the angle
					unless refatom_torsion
						# Should be bound to the atom defining angle
						neighbors = connectivity.bound_atoms(refatom_angle)
						# Discard the ones we do not want
						neighbors.delete(path[i-3])
						neighbors.delete(refatom_dist)
						# Then pick first in the ZM that
						neighbors.each{|ni|
							if zmat.find_atom_index(ni)
								refatom_torsion = ni
								break
							end
						}
					end

					# If there is no reasonable reference for the torsion along the path,
					# define it as improper dihedral.
					unless refatom_torsion
						# Should be bound to the atom defining distance
						neighbors = connectivity.bound_atoms(refatom_dist)

						# Discard the atoms defining the angle
						neighbors.delete(refatom_angle)
						neighbors.delete(atom_index)
						# Then pick first in the ZM that
						neighbors.each{|ni|
							if zmat.find_atom_index(ni)
								refatom_torsion = ni
								break
							end
						}
						raise "Dihedral search failed (improper, possibilities=#{neighbors.size}, linear_angle=#{was_linear})" unless refatom_torsion
					end

					refatom_torsion_zmi = zmat.find_atom_index(refatom_torsion)

					raise "Dihedral search failed" unless refatom_torsion_zmi

					torsion_a = Coordinate.torsion(atom, geometry[refatom_dist], geometry[refatom_angle],  geometry[refatom_torsion]) / Math::PI * 180
					coords << Coord.new(refatom_torsion_zmi, torsion_a, :torsion)
					raise "Torsion refatom not found" unless refatom_torsion_zmi
				end
				zmat.add_line(Line.new(atom.element, coords, atom_index))

			}
		}

		return zmat
	end

	def find_atom_index(geo_index)
		# Translate index in geometry to index in z-matrix
		@lines.each_index{|i|
			if @lines[i].atom_in_geo == geo_index
				return i
			end
		}
		return nil
	end

	def ZMatrix.find_longest_chain(geometry, connectivity)
		# Find longest chain and rank the atoms by distance from one end of the chain
		# Array of atom ranks is returned
		longest_chain_ranks = nil
		max_rank = 0
		geometry.each_index{|i|
			# Start only at terminal atoms
			next unless connectivity.hybridization(i) == 1
			# Do the search
			available_atoms = (0..(geometry.size-1)).to_a
			atom_ranks = []
			current = [i]
			rank = 0
			while available_atoms.size > 0
				# Remove current atoms from available atoms
				available_atoms = available_atoms - current
				# Set rank of current atoms
				current.each{|c|
					atom_ranks[c] = rank
				}
				# Find unused neighbors of current atoms
				neighbors = []
				current.each{|c|
					connectivity.bound_atoms(c).each{|n|
						neighbors << n if available_atoms.include?(n)
					}
				}
				neighbors.uniq!
				# Set current to neighbors,
				# Increase rank for next round
				current = neighbors
				rank += 1
			end
			# Save chain info if its longer than other
			if atom_ranks.max > max_rank
				max_rank = atom_ranks.max
				longest_chain_ranks = atom_ranks
			end
		}
		return longest_chain_ranks
	end

	def ZMatrix.paths_in_tree(geometry, connectivity, graph_tree_ranks, paths = [])
		# Build all paths, starting from the longest
		# If some paths are provided, add to them
		all_paths = paths
		available_atoms = (0..(geometry.size-1)).to_a

		# Remove atoms without ranks
		graph_tree_ranks.each_index{|i|
			if graph_tree_ranks[i].nil?
				available_atoms.delete(i)
			end
		}

		# Remove atoms already in the paths from available
		paths.each{|path|
			available_atoms = available_atoms - path
		}

		while available_atoms.size > 0
			# Pick the available atom with highest rank
			maxrank = -1
			start = nil
			available_atoms.each{|i|
				if graph_tree_ranks[i] > maxrank
					maxrank = graph_tree_ranks[i]
					start = i
				end
			}
			path = []
			atom_i = start
			# Iterate path down to atom ranked 0
			while maxrank >= 0
				# Save atom to path
				path << atom_i
				# Pick a neighbor with rank maxrank - 1
				neighbors = connectivity.bound_atoms(atom_i)
				nextindex = nil
				neighbors.each{|i|
					if graph_tree_ranks[i] == maxrank - 1
						nextindex = i
					end
				}
				atom_i = nextindex
				maxrank -= 1
			end
			available_atoms = available_atoms - path
			all_paths << path.reverse
		end
		return all_paths
	end

	def ZMatrix.outermost_terminal_atom(geometry, connectivity)
		# Working geometry
		geo = geometry
		# Find a center
		center = geo.center
		# Find atom most distant from the center
		outermost_atom = nil
		dist = -1.0
		geo.each_with_index{|atom, i|
			next unless connectivity.hybridization(i) == 1
			d = atom.distance(center)
			if d > dist
				dist = d
				outermost_atom = i
			end
		}
		# Return index of the atom
		return outermost_atom

	end

	def ZMatrix.flood_ranks_add(graph_tree_ranks, geometry, connectivity, subset_list = nil)
		# Builds graph_tree_ranks (into a predefined array which may already contain some rank)
		# starting at atom 'start', possibly operating only on a subset of the geometry
		if subset_list
			available_atoms = subset_list.dup
		else
			available_atoms = (0..(geometry.size-1)).to_a
		end

		# Remove atoms with rank from available
		ranked = []
		available_atoms.each{|i|
			ranked << i if graph_tree_ranks[i]
		}
		available_atoms -= ranked

		# Rank atoms
		while available_atoms.size > 0
			# Find another available atom next to lowest ranking atom
			lowestrank = 9.9e99
			lowranking = []
			available_atoms.each{|i|
				neighbors = connectivity.bound_atoms(i)
				neighbors &= subset_list if subset_list
				neighbors.each{|j|
					if graph_tree_ranks[j]
						if graph_tree_ranks[j] < lowestrank
							lowestrank = graph_tree_ranks[j]
							lowranking = [j]
						elsif graph_tree_ranks[j] == lowestrank
							lowranking << j
						end
					end
				}
			}
			# Grow by 1
			rank = lowestrank + 1
			lowranking.each{|i|
				neighbors = connectivity.bound_atoms(i)
				neighbors &= subset_list if subset_list
				neighbors.each{|j|
					unless graph_tree_ranks[j]
						graph_tree_ranks[j] = rank
						available_atoms.delete(j)
					end
				}
			}
		end
	end

	def ZMatrix.flood_ranks(graph_tree_ranks, geometry, connectivity, start, subset_list = nil)
		# Builds graph_tree_ranks (into a predefined array which may already contain some rank)
		# starting at atom 'start', possibly operating only on a subset of the geometry
		if subset_list
			available_atoms = subset_list.dup
		else
			available_atoms = (0..(geometry.size-1)).to_a
		end
	
		# Rank atoms
		rank = 0
		active = [start]
		while available_atoms.size > 0
			# Rank active atoms and remove them from available atoms
			active.each{|i| 
				graph_tree_ranks[i] = rank
				available_atoms.delete(i)
			}
			# Find all adjacent
			adjacent = []
			active.each{|i|
				neighbors = connectivity.bound_atoms(i) & available_atoms
				adjacent = adjacent | neighbors
			}
			# Move on
			if adjacent.size > 0
				# Expand
				rank += 1
				active = adjacent
			else
				# If none adjacent found, find another available atom next to lowest ranking atom
				lowestrank = 9.9e99
				lowranking = []
				available_atoms.each{|i|
					neighbors = connectivity.bound_atoms(i)
					neighbors &= subset_list if subset_list
					neighbors.each{|j|
						if graph_tree_ranks[j]
							if graph_tree_ranks[j] < lowestrank
								lowestrank = graph_tree_ranks[j]
								lowranking = [j]
							elsif graph_tree_ranks[j] == lowestrank
								lowranking << j
							end
						end
					}
				}
				active = lowranking
				rank = lowestrank + 1
			end
		end
	end
end

