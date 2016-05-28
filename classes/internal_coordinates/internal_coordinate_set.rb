require "classes/internal_coordinates/internal_coordinate.rb"
require "classes/internal_coordinates/jacobian.rb"
require "set"

class InternalCoordinateSet < Array
	# Container for set of internal coordinates
	
	attr_reader :atomlist

	#=======================================================================
	# Constructors
	#=======================================================================
	
	def initialize(geometry, settings = nil)
		super()
		@geometry = geometry
		@atomlist = geometry.atomlist_from_selection(settings[:optimize_region]) if settings and settings[:optimize_region] and settings[:optimize_region] != "%all()"
	end
	
	def InternalCoordinateSet.from_z_matrix(z_matrix, geometry)
		# Traverse Z-matrix, create InternalCoordinateI objects from
		# coordinates defined there
		
		ics = InternalCoordinateSet.new(geometry)
		
		z_matrix.lines.each_index{|i|
			line = z_matrix.lines[i]

			if line.size >= 1
				c = InternalCoordinateI.new(:distance, geometry, line.distance.ref_atom, i)
				c.name = line.distance.name if line.distance.name
				ics << c
			end
			
			if line.size >= 2
				c = InternalCoordinateI.new(:angle, geometry, line.angle.ref_atom, line.distance.ref_atom, i)
				c.name = line.angle.name if line.angle.name
				ics << c
			end

			if line.size >= 3
				c = InternalCoordinateI.new(:torsion, geometry, line.torsion.ref_atom, line.angle.ref_atom, line.distance.ref_atom, i)
				c.name = line.torsion.name if line.torsion.name
				ics << c
			end
		}

		return ics
	end
	
	def InternalCoordinateSet.redundant(geometry, settings)
		ics = InternalCoordinateSet.new(geometry, settings)
		
		# Get covalent connectivity	
		conn = Connectivity.new(geometry)

		# Covalent bonds
		bonds = []
		count = 0
		conn.each_bond{|i,j|
			if (not ics.atomlist) or ics.atomlist.include?(i) or ics.atomlist.include?(j)
				ic = InternalCoordinateI.new(:distance, geometry, i,j)
				ic.name = :bond
				bonds << ic
				count += 1
			end
		}
		Cuby::log.puts_debug("Generating redundant internal coordinates:")
		Cuby::log.puts_debug("   #{count} bonds")
		
		# Add closest intermolecular contacts
		
		# Connect all molecules with van der Waals contact
		count = 0
		molecules = conn.molecules

		if molecules.size > 1
			molecules.each_index{|i|
				i.times{|j|
					# Find closest contact between molecules i and j
					distance, ai, aj = geometry.closest_contact(molecules[i], molecules[j])
					rvdw = PeriodicTable.vdw_radius(geometry[ai].element) + PeriodicTable.vdw_radius(geometry[aj].element)
					if distance < rvdw * 1.2
						if (not ics.atomlist) or ics.atomlist.include?(ai) or ics.atomlist.include?(aj)
							ic = InternalCoordinateI.new(:distance, geometry, ai,aj)
							ic.name = :intermolecular_vdw
							bonds << ic 
							# Add it as a bond to the connectivity,
							# so that adjacent angles and torsions are defined
							count += 1
						end
						conn.add_bond(ai,aj) 
					end
				}
			}

		end
		Cuby::log.puts_debug("   #{count} intermolecular_vdw")

		# Check if there are still isolated molecules
		# connect them to their neighbors
		#-------------------------------------------------- 
		count = 0
		molecules = conn.molecules
		if molecules.size > 1
			Cuby::warning("The molecules can not be connected by links at vdW distance")
		end
		while molecules.size > 1
			dist, mol, ai, aj = geometry.find_closest_molecule(molecules,0)
			if (not ics.atomlist) or ics.atomlist.include?(ai) or ics.atomlist.include?(aj)
				ic = InternalCoordinateI.new(:distance, geometry, ai,aj)
				ic.name = :intermolecular_long
				bonds << ic
				# Add it as a bond to the connectivity,
				# so that adjacent angles and torsions are defined
				count += 1
			end
			conn.add_bond(ai,aj) 
			molecules = conn.molecules
		end
		Cuby::log.puts_debug("   #{count} intramolecular contacts at large distance")


		
		# Add all other close contacts
		#--------------------------------------------------
		# if settings[:redcoord_extra_contacts]
		# 	contact_tol = 1.2
		# 	conn.each_nonbonded_above_1_4{|i,j|
		# 		ai = geometry[i]
		# 		aj = geometry[j]
		# 		if ai.distance(aj) < contact_tol * (PeriodicTable.vdw_radius(ai.element) + PeriodicTable.vdw_radius(aj.element))
		# 			if (not ics.atomlist) or ics.atomlist.include?(i) || ics.atomlist.include?(j)
		# 				bonds << InternalCoordinateI.new(:distance, geometry, i,j)
		# 				conn.add_bond(i,j) # Add the contact to connectivity as well
		# 				Cuby::log.puts_debug("Adding extra distance coordinate between atoms #{i+1} and #{j+1}")
		# 			end
		# 		end
		# 	}
		# end
		#-------------------------------------------------- 


		# When frozen region is used, add bonds adjacent to it so that all bonds and angles can be defined
		# i-|-j-k-l
		#   ^-border
		# add j-k, k-l
		outerbonds = []
		if ics.atomlist
			bonds.each{|bond|
				i,j = bond.indexes
				next if  ics.atomlist.include?(i) && ics.atomlist.include?(j)
				i,j = j,i if ics.atomlist.include?(j) && !ics.atomlist.include?(i)
				# Now, i is inside and j is outside
				# get all k's
				k_list = []
				conn.bound_atoms(j).each{|k|
					next if ics.atomlist.include?(k) # k must be outside, so it can not be i
					k_list << k 
					outerbonds << [j,k].sort # Sorted, so that duplicates can be eliminated easily
				}
				# Get all l's
				k_list.each{|k|
					conn.bound_atoms(k).each{|l|
						next if ics.atomlist.include?(l)
						next if l == j
						outerbonds << [k,l].sort # Sorted, so that duplicates can be eliminated easily
					}
				}

			}
			outerbonds.uniq! # Eliminate duplicates
			# Map it to internal coordinates
			outerbonds.map!{|b| 
				ic = InternalCoordinateI.new(:distance, geometry, b[0], b[1])
				ic.name = :bond_outer
				ic
			}
		end

		outerbonds.each{|b| bonds << b }
		

		# New angles, from bonds
		angles = []
		outerangles = []
		bonds.each_index{|bi1|
			bond1 = bonds[bi1]
			bi1.times{|bi2|
				bond2 = bonds[bi2]
				a,b = bond1.indexes
				c,d = bond2.indexes
				angle = nil
				if a == c
					angle = InternalCoordinateI.new(:angle, geometry, b,a,d)
				elsif a == d
					angle = InternalCoordinateI.new(:angle, geometry, b,a,c)
				elsif b == c
					angle = InternalCoordinateI.new(:angle, geometry, a,b,d)
				elsif b == d
					angle = InternalCoordinateI.new(:angle, geometry, a,b,c)
				end
				
				if angle
					if bond1.name == :bond_outer && bond2.name == :bond_outer
						angle.name = :angle_outer
						outerangles << angle
					end

					if bond1.name.to_s =~ /intermolecular/ || bond2.name.to_s =~ /intermolecular/
						angle.name = :intermolecular
					end

					angles << angle
				end
			}
		}
		Cuby::log.puts_debug("   #{angles.size} angles")
		
		# Torsions - new code, from angles
		torsions = []
		angles.each_index{|ia1|
			angle1 = angles[ia1]
			i,j,k = angle1.indexes
			ia1.times{|ia2|
				angle2 = angles[ia2]
				a,b,c = angle2.indexes

				torsion = nil
				if    j == a && k == b && i != c
					torsion = InternalCoordinateI.new(:torsion, geometry, i,j,k,c)
				elsif j == c && k == b && i != a
					torsion = InternalCoordinateI.new(:torsion, geometry, i,j,k,a)
				elsif j == c && i == b && a != k
					torsion = InternalCoordinateI.new(:torsion, geometry, a,i,j,k)
				elsif j == a && i == b && c != k
					torsion = InternalCoordinateI.new(:torsion, geometry, c,i,j,k)
				end

				if torsion
					if angle1.name == :intermolecular || angle2.name == :intermolecular
						torsion.name = :intermolecular
					end
					torsions << torsion
				end
			}
		}
		Cuby::log.puts_debug("   #{torsions.size} torsions")

		#!# No torsion found
		#!# Improper dihedral should be used
		Cuby::error("No torsion found") if geometry.size >= 4 && torsions.size == 0


		# Connect nonbonded residues in PDB files
		if geometry.is_pdb? && settings[:development][:interfragment_coordinates]
			residues = geometry.residues_as_atomlists
			residues.each_index{|i|
				i.times{|j|
					distance, ai, bi = geometry.closest_contact(residues[i], residues[j])
					next if conn.is_bond?(ai, bi)
					next if distance > 1.2 * (PeriodicTable.vdw_radius(geometry[ai].element) + PeriodicTable.vdw_radius(geometry[bi].element))
					conn.add_bond(ai, bi)
					if (not ics.atomlist) or ics.atomlist.include?(ai) or ics.atomlist.include?(bi)
						ic = InternalCoordinateI.new(:distance, geometry, ai,bi)
						ic.name = :interfragment_vdw
						bonds << ic
					end
				}
			}
		end

		# Write PDB with connectivity
		geometry.write_pdb(:file => "connectivity.pdb", :connectivity => conn) if Cuby::log.logs[0].verbosity == :debug
		
		# Populate the internal coordinate set
		bonds = bonds - outerbonds
		angles = angles - outerangles
		bonds.each{|b| ics << b}
		angles.each{|a| ics << a}
		torsions.each{|t| ics << t}

		# SAMPLE CODE
		# three cartesians for non-frozen molecules, possible way?
		#~ if @settings[:optimize_region].nil? or @settings[:optimize_region] == "%all()" #!# not solving the case of small number of frozen coords
			#~ i,j,k = conn.three_cart_to_hold
			#~ ics << InternalCoordinateI.new(:cartesian, geometry, i)
			#~ ics << InternalCoordinateI.new(:cartesian, geometry, j)
			#~ ics << InternalCoordinateI.new(:cartesian, geometry, k)
		#~ end
		
		return ics
	end
	
	#=======================================================================
	# Calculation of matrices used in conversions
	#=======================================================================
	
	def csr_b
		update_b unless @b_new
		return @csr_b
	end
	
	def csr_bt
		update_b unless @b_new
		return @csr_bt
	end
	
	def csr_btb
		update_b unless @b_new
		return @csr_btb
	end
	
	# updates B-matrices, exploits their structure if they exist
	def update_b
		b_init if @csr_b.nil?
		
		each_index{|row|
			q = at(row)
			cartesians = q.atoms
			q.size.times {|coord_i|
				if (not atomlist) or atomlist.include?(q.indexes[coord_i])
					j = Jacobian.dq_dx(q.type, cartesians, coord_i)
					3.times{|k|
						column = atomlist ? ((atomlist.index(q.indexes[coord_i]))*3+k) : (q.indexes[coord_i]*3+k)
						@csr_b[row, column] = @csr_bt[column, row] = j[k]
					}
				end
			}
		}
		
		@csr_btb.m.times{|i|
			@csr_btb.row_each_nonzero_with_index(i){|x,j|
				@csr_btb[i,j] = @csr_btb[j,i] = @csr_bt.row_dot_row(i,j) if i <= j
		}}
		
		@b_new = true
	end
	
	# initializes structures of B-matrices
	def b_init
		# first construct sparse B-matrix
		b_size = atomlist ? (atomlist.size * 3) : (@geometry.size * 3)
		sparse_b = SparseMatrix.new(self.size, b_size)
		each_index {|row|
			q = at(row)
			cartesians = q.atoms
			q.size.times {|coord_i|
				if (not atomlist) or atomlist.include?(q.indexes[coord_i])
					3.times{|k|
						column = atomlist ? ((atomlist.index(q.indexes[coord_i]))*3+k) : (q.indexes[coord_i]*3+k)
						sparse_b[row, column] = 1.0
					}
				end
			}
		}
		
		# then transform its B, B^T and B^T B variants to CSR format
		@csr_b = SparseMatrixCSX.from_sparse_matrix(:row, sparse_b)
		# tricky transposition
		@csr_bt = SparseMatrixCSX.from_sparse_matrix(:col, sparse_b)
		@csr_bt.change_type!(:row)
		
		sparse_btb = SparseMatrix.new(sparse_b.n)   # empty square sparse matrix
		# diagonal members
		sparse_btb.n.times do |i|
			sparse_btb[i,i] = @csr_bt.row_dot_row(i,i)
		end
		# other membres, here it is O(n^2)
		sparse_btb.n.times do |i|
			for j in i+1 .. sparse_btb.n-1 do
				rd_ij = @csr_bt.row_dot_row(i,j)
				sparse_btb[i,j] = sparse_btb[j,i] = rd_ij unless rd_ij == 0.0
			end
		end
		@csr_btb = SparseMatrixCSX.from_sparse_matrix(:row, sparse_btb)
	end
	
	def geometry_changed
		@b_new = false
		@b_ps_inv = nil
		@p_matrix = nil
		@k_matrix = nil
		return nil
	end
	
	def b
		return self.csr_b.to_matrix
	end
	
	def b_ps_inv
		@b_ps_inv = self.csr_b.to_matrix.pseudoinverse if @b_ps_inv.nil?
		return @b_ps_inv
	end
	
	#=====================================================================
	
	def p_matrix_calc
		@p_matrix = self.b * self.b_ps_inv
	end
	
	def p_matrix
		p_matrix_calc if @p_matrix.nil?
		return @p_matrix
	end
	
	#!# nebude fungovat pro atomlist ani pro nonempty redcoord_freeze_class !!!
	def k_matrix_calc(gradient)
		g_int = self.b_ps_inv.transpose * gradient.to_vector
		k_size = atomlist ? (atomlist.size * 3) : (@geometry.size * 3)
		@k_matrix = Matrix.zero(k_size)
		
		each_index{|i|
			q = at(i)
			cartesians = q.atoms
			q.size.times{|coord_i|
				q.size.times {|coord_j|
					j2 = Jacobian.dq2_dxdy(q.type, cartesians, coord_i, coord_j)
					3.times{|jj|
						3.times{|kk|
							@k_matrix[q.indexes[coord_i]*3+jj, q.indexes[coord_j]*3+kk] += g_int[i] * j2[jj,kk] 
		}}}}}
	end
	
	def k_matrix(gradient)
		k_matrix_calc(gradient) if @k_matrix.nil?
		return @k_matrix
	end
	
	#=======================================================================
	# Misc
	#=======================================================================
	
	def to_s
		s = ""
		each_index{|i|
			s += sprintf("%9s: %3.5f\n", self[i].type.to_s, self[i].value)
		}
		return s
	end
	
	def to_vector
		v = Vector.of_size(self.size)
		each_index{|i|
			v[i] = self[i].value
		}
		return v
	end

	def each_with_index
		each_index{|i|
			yield self[i], i
		}
	end

	def count_types
		counts = {}
		each{|ic|
			if counts[ic.type]
				counts[ic.type] += 1
			else
				counts[ic.type] = 1
			end
		}
		return counts
	end

	def count_type_name
		counts = {}
		each{|ic|
			typename = ic.type.to_s
		       	typename += " (" + ic.name.to_s + ")" if ic.name
			if counts[typename]
				counts[typename] += 1
			else
				counts[typename] = 1
			end
		}
		return counts
	end
	
end
