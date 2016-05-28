class CappedFragments

	class CapAtom
		attr_reader :parent_geometry
		attr_reader :atom_from_i
		attr_reader :atom_orig_i
		attr_reader :dist_ratio

		attr_reader :link_atom

		def initialize(parent_geometry, atom_from_i, atom_orig_i, dist_ratio)
			@parent_geometry = parent_geometry
			@atom_from_i = atom_from_i
			@atom_orig_i = atom_orig_i

			# Original bond scaling ratio, calculated from covalent radii
			if dist_ratio.strip.downcase == "auto"
				@dist_ratio = (PeriodicTable.covalent_radius(:H) + PeriodicTable.covalent_radius(@parent_geometry[@atom_from_i].element)) / 
					(PeriodicTable.covalent_radius(@parent_geometry[@atom_from_i].element) + PeriodicTable.covalent_radius(@parent_geometry[@atom_orig_i].element))
			else
				@dist_ratio = dist_ratio.to_f
			end

			# Create the link atom
			@link_atom = Atom.new(:H)
			@link_atom.properties[:pdb_atom_name] = "HL"
		end

		def to_s
			return "Cap, H, g= #{@dist_ratio}"
		end

		def update_position!
			atom_from = @parent_geometry[@atom_from_i]
			atom_to = @parent_geometry[@atom_orig_i]
			@link_atom.set_coord((atom_to - atom_from)* @dist_ratio + atom_from)
		end
	end

	class Fragment
		attr_reader :order
		attr_reader :atomlist_in_parent
		attr_reader :parent_geometry
		attr_reader :id

		attr_reader :caps
		attr_reader :geometry

		attr_accessor :calculation

		ORDER_NAMES = {
			1 => :monomer,
			2 => :dimer,
			3 => :trimer,
			4 => :tetramer
		}

		def initialize(order, atomlist, parent_geometry, id, caps)
			@order = order
			@atomlist_in_parent = atomlist
			@parent_geometry = parent_geometry
			@id = id.to_s

			# Create geometry
			@geometry = @parent_geometry.geometry_from_list(@atomlist_in_parent)

			# Caps
			@caps = caps
			# Add caps to geometry
			@caps.each{|cap|
				@geometry << cap.link_atom
			}
			
		end

		def Fragment.merge(frag1, frag2)
			# Merge atomlists
			atomlist = frag1.atomlist_in_parent | frag2.atomlist_in_parent

			# Merge caps - eliminate these within merged dimer
			merged_caps = []
			(frag1.caps | frag2.caps).each{|cap|
				merged_caps << cap unless atomlist.include?(cap.atom_orig_i)
			}

			# Merge the monomers
			merged = Fragment.new(
				frag1.order + frag2.order,
				atomlist,
				frag1.parent_geometry,
				"#{frag1.id}-#{frag2.id}",
				merged_caps
			)

			return merged
		end

		def to_s
			return sprintf("%10s: %10s, %3d atoms, %3d caps", @id, ORDER_NAMES[@order].to_s, @atomlist_in_parent.size, @caps.size)
		end

		def update_caps_position!
			@caps.each{|cap| cap.update_position!}
		end
		
	end

	attr_reader :parent_geometry

	attr_reader :monomers
	attr_reader :dimers

	def initialize(parent_geometry, settings)
		@parent_geometry = parent_geometry

		# Do the fragmentation
		# yields geometries of fragments with data on their interconnection
		case settings[:fragmentation_mode]
		when :cut_bonds
			@monomers = build_fragments_cut_bonds(settings)
		end


		# Set position of link atoms
		update_link_atoms

		# Build dimers
		@dimers = build_dimers

		# Save geometries to PFB file
		save_fragments_pdb(settings[:fragmentation_monomer_file]) unless settings[:fragmentation_monomer_file] == ""
		save_dimers_pdb(settings[:fragmentation_dimer_file]) unless settings[:fragmentation_dimer_file] == ""
	end

	def all_fragments
		return @monomers | @dimers.flatten
	end

	#=======================================================================
	# Methods
	#=======================================================================
	
	def build_fragments_cut_bonds(settings)
		# Build fragments by cutting bonds specified in settings

		# Shortcut:
		geo = @parent_geometry

		# Build connectivity
		connectivity = Connectivity.new(geo)

		# Expand the selections
		# The selections can be generic, make list of real bonds matching all combinations
		expanded_bond_list = {}
		settings[:fragmentation_cut_bonds].each{|line|
			Cuby::error("The entries of the list of bonds to be cut must be in a form \"selection;selection;distance_ratio\"") unless line =~ /[^;]+;[^;]+;[^;]+/
			sel1, sel2, ratio = line.strip.split(/\s*;\s*/)
			list1 = geo.atomlist_from_selection(sel1)
			list2 = geo.atomlist_from_selection(sel2)
			list1.each{|at1|
				list2.each{|at2|
					if connectivity.is_bond?(at1, at2)
						# Save bond, in both directions
						expanded_bond_list[at1] = [at2, ratio]
						expanded_bond_list[at2] = [at1, ratio]
					end
				}
			}
		}
		fragments = connectivity.molecules
		Cuby::log.puts_v(:normal, "Number of non-covalent fragments: #{fragments.size}")
		Cuby::log.puts_v(:normal, "Number of bonds cut: #{expanded_bond_list.size.to_f/2}")

		# Make cuts
		cut_bonds = {}
		expanded_bond_list.each_pair{|at1, item|
			at2, ratio = item
			# Erase from connectivity
			connectivity[at1,at2] = 0
		}

		# Get fragments geometries
		fragments = connectivity.molecules # Yields atomlist
		Cuby::log.puts_v(:normal, "Number of resulting fragments: #{fragments.size}")

		monomers = []
		fragments.each_index{|i|
			fragment = fragments[i] # Atomlist

			# Find caps
			caps = []
			fragment.each{|at_i|
				if expanded_bond_list[at_i]
					# Cut bond starts from here
					cap = CapAtom.new(@parent_geometry, at_i, expanded_bond_list[at_i][0], expanded_bond_list[at_i][1])
					caps << cap
				end
			}

			# Create monomer
			monomers[i] = Fragment.new(1, fragment, @parent_geometry, i, caps)
		} 

		# Print fragments
		#--------------------------------------------------
		# puts "Monomers:"
		# monomers.each{|mono| puts mono}
		#-------------------------------------------------- 

		return monomers
	end

	def update_link_atoms
		# Set positions of link atoms
		@monomers.each{|mono|
			mono.update_caps_position!
		}

	end

	def save_fragments_pdb(filename)
		# Save fragments as residues to a PDB file
		geo = Geometry.new
		@monomers.each_index{|i|
			f2 = @monomers[i].geometry.deep_copy
			f2.each{|atom|
				atom.properties[:pdb_res_no] = i
				atom.properties[:pdb_res_name] = "f%02d" % (i+1)
				geo << atom
			}
			f2.last.properties[:terafter] = true
		}
		geo.write_pdb(:file => filename, :long => filename =~ /\.lpdb$/)
	end

	def build_dimers # => array of geometries
		# Builds all dimers
		dimers = []
		@monomers.each_index{|i|
			dimers[i] = []
			frag_i = @monomers[i]
			i.times{|j|
				frag_j = @monomers[j]
				dimers[i][j] = Fragment.merge(frag_i, frag_j)
			}
		}

		# Print
		#--------------------------------------------------
		# puts "Dimers:"
		# dimers.flatten.each{|dimer| puts dimer}
		#-------------------------------------------------- 

		return dimers
	end

	def save_dimers_pdb(filename)
		# Save dimers as residues to a PDB file
		geo = Geometry.new
		count = 1
		@dimers.each_index{|i|
			i.times{|j|
				f2 = @dimers[i][j].geometry.deep_copy
				f2.each{|atom|
					atom.properties[:pdb_res_no] = count
					atom.properties[:pdb_res_name] = "d%02d" % count
					geo << atom
				}
				f2.last.properties[:terafter] = true
				count += 1
			}
		}
		geo.write_pdb(:file => filename, :long => filename =~ /\.lpdb$/)
	end
end
