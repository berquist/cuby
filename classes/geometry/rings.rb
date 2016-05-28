################################################################################
#
# Module GeometryRings
#
# Author: Jan Rezac
# Date created: 2010-08-18
# License: Cuby license
# Description: Search for rings in the geometry
# Extends: Geometry
# Status: Works
#
################################################################################

require "classes/algebra/matrix_graph"

class Ring < Array
	# Ring defined as a list of atom indexes
	attr_reader :parent_set

	def initialize(parent_set)
		@parent_set = parent_set
		super()
	end

	def is_aromatic?
		#--------------------------------------------------
		# # not aromatic if not planar
		# if ring.size >= 4
		# 	#!# Does not work, must be done better - get ring center, average plane norm vector over all the bonds...
		# 	a = @geometry.at(ring[0])
		# 	b = @geometry.at(ring[1])
		# 	c = @geometry.at(ring[2])
		# 	plane_n = (a-b).cross_product(c-b)
		# 	(ring.size - 3).times{|i|
		# 		dplane = @geometry.at(i+3).point_plane_distance(plane_n, b)
		# 		puts dplane
		# 		return false if dplane > 0.2
		# 	}
		# end
		#-------------------------------------------------- 
		geometry = @parent_set.geometry
		conn = @parent_set.connectivity
		pi_electrons = 0
		each{|i|
			atom = geometry.at(i)
			bonds = conn.hybridization(i)
			case atom.element
			when :C
				if bonds == 3
					# find the atom outside the ring
					ao = nil
					conn.bound_atoms(i).each{|a| 
						ao = a unless self.include?(a)
					}
					if geometry[ao].element == :O && conn.hybridization(ao) == 1
						# no pi electrons if there is a double bond out of the ring (look for =O)
					else
						pi_electrons += 1
					end
				end
			when :N
				pi_electrons += 1 if bonds == 2
			when :O
				pi_electrons += 2 if bonds == 2
			else
				raise "Aromaticity detection works only for C,N,O"
			end
		}
		return true if (pi_electrons - 2) % 4 == 0
		return false
	end
end

class Rings

	attr_reader :geometry
	attr_reader :connectivity

	def initialize(geometry, size_limit = 9999999)
		#: Search for Smallest Set of Smallest Rings (SSSR). In SSSR, no ring is
		#: superset of a smaller ring.

		# For reference, see MatrixGraph.smallest_set_of_smallest_rings
		
		# Save geometry
		@geometry = geometry
		if @geometry.info[:connectivity]
			@connectivity = @geometry.info[:connectivity]
		else
			@connectivity = Connectivity.new(@geometry)
		end

		# Reduce to core
		#!# use @connectivity rather than recalculating it in the following method
		core = reduce_to_ring_core(@geometry)

		@rings = []
		if core.size > 0
			# Build graph
			graph = connectivity_to_graph(core, 1.0)

			# Get SSSR of the graph
			sssr = graph.smallest_set_of_smallest_rings(size_limit)

			# Translate rings to atom indices in the original system
			sssr.each_index{|i|
				ring = sssr[i]
				atomlist = Ring.new(self)
				ring.each{|v|
					atomlist << @geometry.index(core[v])
				}
				@rings[i] = atomlist
			}
		end

		return nil
	end

	def [] (index)
		return @rings[index]
	end

	def  size
		return @rings.size
	end

	def rings_containing_atom(atom_index)
		#: List of rings containing an atom
		list = []
		@rings.each{|ring|
			list << ring if ring.include?(atom_index)
		}
		return list
	end

	def aromatic_atom?(atom_index)
		rings_containing_atom(atom_index).each{|ring|
			return true if ring.is_aromatic?
		}
		return false
	end

	#=======================================================================
	# Private
	#=======================================================================

	def reduce_to_ring_core(geometry) # => Geometry
		#: Creates new geometry that have all non-cyclic branches removed. The ring search
		#: is performed on this core only.
		work_geo = geometry
		n = work_geo.size
		n_prev = n + 1
		while n != n_prev
			conn = Connectivity.new(work_geo)
			new_geo = Geometry.new
			work_geo.each_with_index{|atom, i|
				new_geo << atom unless conn.hybridization(i) <= 1.0
			}
			work_geo = new_geo
			n_prev = n
			n = work_geo.size
			break if n == 0
		end
		return work_geo
	end

	def connectivity_to_graph(geometry, bond_weight) # => MatrixGraph
		#: Converts connecivity information to graph.
		conn = Connectivity.new(geometry)
		graph = MatrixGraph.empty(geometry.size)
		geometry.size.times{|i|
			i.times{|j|
				graph.add_edge_undirected(i,j, bond_weight) if conn.is_bond?(i,j)
			}
		}
		return graph
	end
end
