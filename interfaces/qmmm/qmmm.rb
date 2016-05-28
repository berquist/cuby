################################################################################
#
# QMMM interface
#
# Author: Jan Rezac
# Date created: 2014-01-23
# License: Cuby4 license
# Description: QM/MM interface
# Status: Works
#
################################################################################

#===============================================================================
# Uses subtractive scheme to calculate interaction across the QM/MM boundary
# Gradient projection from:
# Dapprich S., Komaromi I., Suzie Byun K., Morokuma K., Frisch M.J.,
# J. Mol. Struct. Teochem 1999
#===============================================================================

# ToDo:
# * Hessian (cuby3)
# * Support for external point charges for the whole setup (cuby3)
# * Microiteraions (cuby3)

module InterfaceQmmm
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :composite
	CAPABILITIES = [:energy, :gradient]
	MODIFIER = false
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation_qm,          :required, "Setup for the QM calculation"],
		InputBlock[:calculation_mm,          :required, "Setup for the MM calculation"],
		InputBlock[:calculation_qmregion_mm, :optional, "MM setup applied only to the QM region"],
		InputBlock[:calculation_system_mm,   :optional, "MM setup applied only to the calulation of the whole system"]
	]
	#=======================================================================

	def prepare_interface
		# Checks
		if  @settings[:qmmm_embedding] == :electrostatic && @what.include?(:gradient)
			Cuby::warning("Gradients obtained with qmmm_embedding == electrostatic are not exact,\nuse electrostatic_2way instead.")
		end

		# Auto fragmentation
		case @settings[:qmmm_auto_fragmentation]
		when :peptide_backbone
			auto_fargmentation_peptide_backbone
		end
		
		# Prepare data needed, build geometry of the QM region
		qmmm_fill_lists	
		qmmm_build_qm_region
		qmmm_update_coords

		# Save geometry
		unless @settings[:qmmm_qmregion_file] == ""
			Cuby::log.puts_v(:normal, "Writting QM region geometry to #{@settings[:qmmm_qmregion_file]}")
			begin
				filetype = GeoFile.type_from_filename(@settings[:qmmm_qmregion_file])
			rescue
				filetype = :pdb
			end
			@qm_region.write_file(@settings[:qmmm_qmregion_file], filetype, {:extra_columns => @settings[:pdb_extra_columns]})
		end

		# Exit if only geometry generation was requested
		if @settings[:qmmm_geometry_only]
			if @settings[:prepare_only]
				Cuby::log.puts_debug("qmmm_geometry_only set, skipping prepare")
				return
			else
				Cuby::log.puts_debug("qmmm_geometry_only set, exitting")
				exit 
			end
		end

	 	# Build calculations
		qmmm_build_calculations
	end

	def queue_interface(queue)
		# Update coordinates of the QM region
		qmmm_update_coords

		# Update coordinates of the point charges
		if @settings[:qmmm_embedding] != :mechanical
			update_point_charges(@calculation_qm_region.point_charges)
		end

		@calculations.each{|calculation|
			calculation.send_to_queue(queue)
		}
	end

	def compose_interface
		return qmmm_results
	end

	def cleanup_interface
		@calculations.each{|calculation|
			calculation.cleanup
		}
	end

	#=======================================================================
	# Private methods: preparation
	#=======================================================================
	
	class QmmmLinkAtom
		attr_accessor :qm_index
		attr_accessor :mm_index
		attr_accessor :ratio
		attr_accessor :pdb_name

		def initialize(qm_index, mm_index)
			@qm_index = qm_index
			@mm_index = mm_index
		end
	end
	
	def qmmm_fill_lists
		# List of atoms in QM core (without link atoms)
		@qmmm_core_list = @geometry.atomlist_from_selection(@settings[:qmmm_core])

		# Link atoms
		@qmmm_link_atoms = []
		@settings[:qmmm_cut_bonds].each_index {|i|
			# Check record format
			bond = @settings[:qmmm_cut_bonds][i]
			Cuby::error "Bond specification (item of qmmm_cut_bonds array) should be of type Hash" if  bond.class != Hash
			Cuby::error "Bond information missing in record #{i} of qmmm_cut_bonds" if bond['bond'].class != String
			Cuby::error "Bond information in record #{i} of qmmm_cut_bonds has wrong format" if bond['bond'] !~ /[0-9]+-[0-9]+/
			# Get indexes of atoms in the bond
			atoms = bond['bond'].split('-').collect{|s| s.to_i - 1}
			atoms_in_core = []
			atoms_not_in_core = []
			atoms.each {|a|
				if @qmmm_core_list.include?(a)
					atoms_in_core << a
				else
					atoms_not_in_core << a
				end
			}
			# Check if it is bond across boundary
			unless atoms_in_core.size == 1 && atoms_not_in_core.size == 1
				Cuby::error "Bond in record #{i} of qmmm_cut_bonds is not at QM/MM boundary"
			end
			# Create new link atom record
			linkatom = QmmmLinkAtom.new(atoms_in_core[0], atoms_not_in_core[0])
			@qmmm_link_atoms[i] = linkatom

			# Save link atom bond distance
			if bond['link_ratio'].kind_of?(Numeric)
				linkatom.ratio = bond['link_ratio'].to_f
			elsif bond['link_dist'].kind_of?(Numeric)
				# calculate ratio from distances
				Cuby::warning "Link atom distance set directly using 'link_dist' option, the use of 'link_ratio' is recommended."
				ratio = bond['link_dist'].to_f / @geometry[linkatom.qm_index].distance(@geometry[linkatom.mm_index])
				linkatom.ratio = ratio
			else
				Cuby::error "Link distance or ratio missing in record #{i} of qmmm_cut_bonds"
			end

			# Save link atom type
			if bond['link_type']
				Cuby::error "Link atom type is not a string in record #{i} of qmmm_cut_bonds" if bond['link_type'].class != String
				linkatom.pdb_name = bond['link_type']
			else
				Cuby::warning "Link atom type not specified, using default value 'HL'"
				linkatom.pdb_name = 'HL'
			end

		}
	end

	def qmmm_build_qm_region
		# QM region geometry
		@qm_region = Geometry.new
		# Indexes of the atoms in:
		# whole system - positive
		# link atom list - negative - 1
		@qm_region_i = []
		# Atomic charges on the atoms
		@qm_region_ch = AtomicCharges.new(:forcefield) if @atomic_charges

		# Temporary list of link atom indexes
		link_list = []
		@qmmm_link_atoms.each_index{|i|
			link_list[i] = @qmmm_link_atoms[i].mm_index
		}

		@geometry.each_index{|i|
			if @qmmm_core_list.include?(i)
				# An atom in the QM region
				@qm_region << @geometry[i].deep_copy
				@qm_region_i << i
				@qm_region_ch << @atomic_charges[i] if @atomic_charges
			elsif link_list.include?(i)
				# MM atom from which a link atoms is constructed
				link_rec_index = link_list.index(i)
				# Create new atom
				link = Atom.new(:H)
				# Set its PDB properties
				link.properties[:pdb_res_no] = @geometry[i].properties[:pdb_res_no] if @geometry[i].properties[:pdb_res_no]
				link.properties[:pdb_res_name] = @geometry[i].properties[:pdb_res_name] if  @geometry[i].properties[:pdb_res_name]
				link.properties[:pdb_atom_name] = @qmmm_link_atoms[link_rec_index].pdb_name
				# Save it to the lists
				@qm_region << link
				@qm_region_i << -link_rec_index - 1
				@qm_region_ch << 0.0 if @atomic_charges
			end
		}

		# Rename qm_region residues
		@settings[:qmmm_rename_residues].each{|res|
			selection, name = res.split
			@qm_region.atomlist_from_selection(selection).each{|i|
				@qm_region[i].properties[:pdb_res_name] = name
			}
		}

		@settings[:qmmm_add_ter].split(/\s;\s/).each{|res|
			@qm_region.geometry_from_selection(res).each{|atom|
				atom.properties[:terafter] = true
			}
		}

		return nil
	end

	def qmmm_update_coords
		@qm_region.each_index{|i|
			atom = @qm_region[i]
			gi = @qm_region_i[i]
			if gi >= 0
				# It is a normal atom
				atom.set_coord(@geometry[gi])
			else
				# Link atom
				link_info = @qmmm_link_atoms[-gi-1]
				atom.set_coord(@geometry[link_info.qm_index] + link_info.ratio * (@geometry[link_info.mm_index] - @geometry[link_info.qm_index]))
			end
		}
	end

	def qmmm_build_calculations
		# At least two blocks are expected: calculation_qm and calculation_mm
		unless @settings.has_block?(:calculation_qm)
			Cuby::error "The QM calculation must be defined in a block 'calculation_qm'"
		end
		unless @settings.has_block?(:calculation_mm)
			Cuby::error "The MM calculation must be defined in a block 'calculation_mm'"
		end

		# Blocks :calculation_qmregion_mm and :calculation_system_mm hold
		# the setup of the MM calculations. If not present, they are created
		@settings.new_block(:calculation_qmregion_mm) unless @settings.has_block?(:calculation_qmregion_mm)
		@settings.new_block(:calculation_system_mm) unless @settings.has_block?(:calculation_system_mm)

		@settings.block(:calculation_qmregion_mm).copy_from!(@settings.block(:calculation_mm), :keep)
		@settings.block(:calculation_system_mm).copy_from!(@settings.block(:calculation_mm), :keep)

		# Build the calculations
		@calculation_qm_region = Calculation.new(@name + "_qmregion_qm", @settings.block(:calculation_qm), @qm_region)
		@calculation_mm_region = Calculation.new(@name + "_qmregion_mm", @settings.block(:calculation_qmregion_mm), @qm_region)
		@calculation_mm_system = Calculation.new(@name + "_system_mm", @settings.block(:calculation_system_mm), @geometry)

		# What to calculate
		system_what = @what.dup
		region_what = @what.dup

		# Prepare the MM calculation - needed to get the point charges
		@calculation_mm_system.prepare(system_what)

		# Prepare point charges for electrostatic embedding
		if @settings[:qmmm_embedding] != :mechanical
			qmmm_prepare_point_charges 
			point_charges = qmmm_get_point_charges
			@calculation_qm_region.point_charges = point_charges
			@calculation_mm_region.point_charges = point_charges
		end
		if @settings[:qmmm_embedding] == :electrostatic_2way
			region_what << :point_charges_gradient
		end

		# Prepare the QM region calculations
		@calculation_qm_region.prepare(region_what)
		@calculation_mm_region.prepare(region_what)

		@calculations = [
			@calculation_qm_region,
			@calculation_mm_region,
			@calculation_mm_system
		]
	end

	def qmmm_prepare_point_charges
		if @settings.set?(:atomic_charges_read)
			# Read the atomic charges from a file
			all_charges = AtomicCharges.from_file(@settings[:atomic_charges_read])
		else
			# Try to get atomic charges from the MM calculation
			if @calculation_mm_system.respond_to?(:get_atomic_charges)
				all_charges = @calculation_mm_system.get_atomic_charges
			else
				Cuby::error "Electrostatic embedding requires atomic charges for all atoms.
these can be either provided from a file (keyword 'atomic_charges_read')
or obtained automatically from MM interface that supports this feature."
			end
		end

		# Point charges list: hash atom_index => charge
		@point_charges_list = {}

		# Charges removed as specified in the input
		user_removed = @geometry.atomlist_from_selection(@settings[:qmmm_remove_charges])

		# Charges on atoms replaced by a link atom
		cutbond_list = []
		@qmmm_link_atoms.each_index{|i|
			cutbond_list[i] = @qmmm_link_atoms[i].mm_index
		}

		# Charges around link atoms
		around_links = []
		if @settings[:qmmm_charges_around_links] > 0
			around_links = list_of_atoms_connected_to_link(@settings[:qmmm_charges_around_links], cutbond_list)
		end

		indexes =  (0..(@geometry.size-1)).to_a - @qmmm_core_list - cutbond_list - user_removed - around_links
		indexes.each{|i|
			@point_charges_list[i] = all_charges[i]
		}

		# Debugging: write point charges to a file
		if Cuby::log.logs[0].verbosity == :debug
			chg = @geometry.geometry_from_list(indexes)
			chg.write_xyz(:file => "qmmm_point_charges.xyz")
		end

		return nil
	end

	def list_of_atoms_connected_to_link(number_of_bonds, starting_indexes)
		connectivity = Connectivity.new(@geometry)
		atoms = starting_indexes.dup
		number_of_bonds.times{
			boundary = []
			atoms.each{|i|
				connectivity.row_each_nonzero_index(i){|j|
					boundary |= [j]
				}
			}
			atoms |= boundary
		}
		return atoms
	end

	def qmmm_get_point_charges
		# build embedding charges for current geometry
		point_charges = PointCharges.new
		@point_charges_list.each_pair{|i, charge|
			pch = PointCharge.from_coordinate(@geometry[i], charge)
			pch.from_atom = @geometry[i] if @settings[:qmmm_charges_extra_info]
			point_charges << pch
		}
		return point_charges
	end

	def update_point_charges(point_charges)
		# Update coordinates of the point charges
		chi = 0
		@point_charges_list.each_pair{|i, charge|
			point_charges[chi].set_coord(@geometry[i])
			chi += 1
		}
		return nil
	end

	#=======================================================================
	# Private methods: results
	#=======================================================================
	
	def qmmm_results
		results = Results.new

		results.energy = @calculation_mm_system.results.energy -
			@calculation_mm_region.results.energy +
			@calculation_qm_region.results.energy

		# Components
		results.energy_components[:qmmm_qmregion_qm] = @calculation_qm_region.results.energy
		results.energy_components[:qmmm_qmregion_mm] = @calculation_mm_region.results.energy
		results.energy_components[:qmmm_system_mm] = @calculation_mm_system.results.energy

		#!# Components of components
		
		# Gradient 
		results.gradient = qmmm_gradient if @what.include?(:gradient)

		# Hessian
		results.hessian = qmmm_hessian if @what.include?(:hessian)

		return results
	end

	def qmmm_gradient
		# Start with gradient of the whole system
		gradient = @calculation_mm_system.results.gradient

		# Add the QM "corerction"
		@qm_region.each_index{|i|
			atom = @qm_region[i]
			gi = @qm_region_i[i]
			if gi >= 0
				# It is a normal atom
				gradient[gi].minus!(@calculation_mm_region.results.gradient[i])
				gradient[gi].plus!(@calculation_qm_region.results.gradient[i])
			else
				# Link atom
				link_info = @qmmm_link_atoms[-gi-1]
				g = link_info.ratio
				dgrad = @calculation_qm_region.results.gradient[i] - @calculation_mm_region.results.gradient[i]
				# Correct the QM atom
				gradient[link_info.qm_index].plus!(dgrad * (1.0-g))
				# Correct the MM atom
				gradient[link_info.mm_index].plus!(dgrad * g)
			end
		}

		# Two-way electrostatic embedding: correct gradient on MM atoms
		if @settings[:qmmm_embedding] == :electrostatic_2way
			chi = 0
			@point_charges_list.each_pair{|i, charge|
				gradient[i].minus!(@calculation_mm_region.results.point_charges_gradient[chi])
				gradient[i].plus!(@calculation_qm_region.results.point_charges_gradient[chi])
				chi += 1
			}
		end

		return gradient
	end

	def qmmm_hessian
	end

	#=======================================================================
	# Auto fragmentation, rebuilds input
	#=======================================================================

	def auto_fargmentation_peptide_backbone
		# Check geometry - it must contain the PDB data
		unless @geometry[0].properties[:pdb_atom_name]
			Cuby::error("QMMM auto fragmentation possible only on geometry from a PDB file")
		end

		# Save original selection
		selection = @settings[:qmmm_core]
		Cuby::log.puts "=" * 80
		Cuby::log.puts "Automated QMMM fragmentation of peptide backbone"
		Cuby::log.puts "=" * 80

		# Save table of residue names
		residue_names = {}
		geometry.each{|atom|
			residue_names[atom.properties[:pdb_res_no]] = atom.properties[:pdb_res_name]
		}

		# Convert the selection to a list of residues
		residues = []
		geometry.atomlist_from_selection(selection).each{|i|
			unless residues.include?(geometry[i].properties[:pdb_res_no])
				residues << geometry[i].properties[:pdb_res_no] 
			end
		}
		residues.sort!
		Cuby::log.puts "Residues in the QM region:"
		#Cuby::log.puts "   " + integer_list_to_string(residues, 0, ", ")
		Cuby::log.puts "   " + residues.join(', ')

		# Fill 1AA gaps
		add_res = []
		caps_b = []
		residues.each_index{|i|
			next if i == 0
			if residues[i] == residues[i-1] + 2
				case residue_names[residues[i-1] + 1] 
				when "GLY"
					# Add as it is
					add_res << residues[i-1] + 1
				when "PRO"
					# Add as it is
					add_res << residues[i-1] + 1
				else
					### Add backbone fragment
					#add_res << residues[i-1] + 1
					caps_b << residues[i-1] + 1
				end
			end
		}
		if add_res.size > 0
			Cuby::log.puts "Adding whole residues to fix the holes:"
			Cuby::log.puts "   " + add_res.join(', ')
		end
		if caps_b.size > 0
			Cuby::log.puts "Adding backbone fragments to fix the holes:"
			Cuby::log.puts "   " + caps_b.join(', ')
		end

		residues = (residues | add_res).sort

		Cuby::log.puts "Used residues:"
		Cuby::log.puts "   " + residues.join(', ')

		# List of residues to be transformed to caps
		caps_n = []
		caps_c = []

		res_count = geometry.pdb_count_residues(true)
		res_count.times{|ii|
			i = ii + 1
			next if i == 1
			next if i == res_count
			if (residues | caps_b).include?(i) && ! (residues | caps_b).include?(i-1)
				# Check for TER before the residue, if found, no cap needed
				unless geometry.geometry_from_selection(":#{i-1}").last.properties[:terafter]
					caps_n << i-1
				end
			end
			if (residues | caps_b).include?(i) && ! (residues | caps_b).include?(i+1)
				# Check for TER at the end of the selected residue
				unless geometry.geometry_from_selection(":#{i}").last.properties[:terafter]
					caps_c << i+1
				end
			end
		}

		Cuby::log.puts "Residues changed to N caps:"
		Cuby::log.puts "   " + caps_n.join(', ')
		Cuby::log.puts "Residues changed to C caps:"
		Cuby::log.puts "   " + caps_c.join(', ')

		# Process caps
		cap_orig_atoms = []
		qmmm_cut_bonds = []
		qmmm_rename_residues = []
		qmmm_batoms = []

		# Process backbone caps
		caps_b.each{|cap|
			# Add atoms to the list
			atoms = {}
			["C","O","CA","CB","HA","N","H"].each{|atomname|
				atom = geometry.atomlist_from_selection("%atomname(#{atomname})&:#{cap}")[0]
				unless atom
					Cuby::error("Can not change residue #{cap} to BCB, atom #{atomname} not found")
				end
				cap_orig_atoms << atom unless atomname == "CB"
				atoms[atomname] = atom
			}
			# Residue renaming
			qmmm_rename_residues << "  - \":#{cap} BCB\""
			# Bond breaking
			# ### Change to optimized value
			qmmm_cut_bonds << "  - { bond: #{atoms["CA"]+1}-#{atoms["CB"]+1}, link_ratio: 0.723, link_type: HB }"
			qmmm_batoms << atoms["CB"]
		}

		# Process N caps
		caps_n.each{|cap|
			#resgeo = geometry.geometry_from_selection(":#{cap}")
			#resgeo.write_pdb
			atom_c  = geometry.atomlist_from_selection("%atomname(C)&:#{cap}")[0]
			atom_o  = geometry.atomlist_from_selection("%atomname(O)&:#{cap}")[0]
			atom_ca = geometry.atomlist_from_selection("%atomname(CA)&:#{cap}")[0]
			Cuby::error "Residue #{cap} can not be changed into a cap, atom C not found" unless atom_c
			Cuby::error "Residue #{cap} can not be changed into a cap, atom O not found" unless atom_o
			Cuby::error "Residue #{cap} can not be changed into a cap, atom CA not found" unless atom_ca
			# Add atoms to the list
			cap_orig_atoms << atom_c
			cap_orig_atoms << atom_o
			# Residue renaming
			qmmm_rename_residues << "  - \":#{cap} BNC\""
			# Bond breaking
			# PM6 distance: 1.1008
			# AMBER CT-C:  1.522
			# ratio in PM6/AMBER FF03 = 0.662
			qmmm_cut_bonds << "  - { bond: #{atom_ca+1}-#{atom_c+1}, link_ratio: 0.723, link_type: HA }"
			qmmm_batoms << atom_ca
		}

		# Process C caps
		caps_c.each{|cap|
			#resgeo = geometry.geometry_from_selection(":#{cap}")
			#resgeo.write_pdb
			atom_n  = geometry.atomlist_from_selection("%atomname(N)&:#{cap}")[0]
			atom_h  = geometry.atomlist_from_selection("%atomname(H)&:#{cap}")[0]
			atom_ca = geometry.atomlist_from_selection("%atomname(CA)&:#{cap}")[0]
			atom_ha = geometry.atomlist_from_selection("%atomname(HA)&:#{cap}")[0]
			atom_cb = geometry.atomlist_from_selection("%atomname(CB)&:#{cap}")[0]
			atom_c  = geometry.atomlist_from_selection("%atomname(C)&:#{cap}")[0]

			resname = geometry[geometry.atomlist_from_selection(":#{cap}")[0]].properties[:pdb_res_name]
			case resname 
			when 'PRO'
				# Proline, does not work
				Cuby::error "Residue #{cap} (#{resname}) can not be changed into a cap, add this residue to QM part selection"
			when 'GLY'
				# Glycine, BCG cap
				atom_2ha = geometry.atomlist_from_selection("(%atomname(2HA)|%atomname(HA2))&:#{cap}")[0]
				atom_3ha = geometry.atomlist_from_selection("(%atomname(3HA)|%atomname(HA3))&:#{cap}")[0]
				# Add atoms to the list
				cap_orig_atoms << atom_n
				cap_orig_atoms << atom_h
				cap_orig_atoms << atom_ca
				cap_orig_atoms << atom_2ha
				cap_orig_atoms << atom_3ha
				# Residue renaming
				qmmm_rename_residues << "  - \":#{cap} BCG\""
				# PM6 distance: 1.103
				# AMBER CT-C: 1.522
				# ratio in PM6/AMBER FF03 = 0.725
				qmmm_cut_bonds << "  - { bond: #{atom_ca+1}-#{atom_c+1}, link_ratio: 0.725, link_type: H2 }"
				qmmm_batoms << atom_c
			else
				# Other residues
				Cuby::error "Residue #{cap} (#{resname}) can not be changed into a cap, atom CB not found" unless atom_cb
				Cuby::error "Residue #{cap} (#{resname}) can not be changed into a cap, atom H not found" unless atom_h
				# Add atoms to the list
				cap_orig_atoms << atom_n
				cap_orig_atoms << atom_h
				cap_orig_atoms << atom_ca
				cap_orig_atoms << atom_ha
				# Residue renaming
				qmmm_rename_residues << "  - \":#{cap} BCC\""
				# Bond breaking
				# PM6 distance: 1.103
				# AMBER CT-CT: 1.526
				# ratio in PM6/AMBER FF03 = 0.723
				qmmm_cut_bonds << "  - { bond: #{atom_ca+1}-#{atom_cb+1}, link_ratio: 0.723, link_type: HB }"
				# PM6 distance: 1.103
				# AMBER CT-C: 1.522
				# ratio in PM6/AMBER FF03 = 0.725
				qmmm_cut_bonds << "  - { bond: #{atom_ca+1}-#{atom_c+1}, link_ratio: 0.725, link_type: H2 }"
				qmmm_batoms << atom_c
				qmmm_batoms << atom_cb
			end
		}

		# Copy this setup into settings
		@settings[:qmmm_core] = ":#{integer_list_to_string(residues, 0, ",")}|@#{integer_list_to_string(cap_orig_atoms,1)}"
		@settings[:qmmm_cut_bonds] = []
		qmmm_cut_bonds.each{|s|
			@settings[:qmmm_cut_bonds] << YAML.load(s.gsub(/^ *- */,''))
		}
		@settings[:qmmm_rename_residues] = []
		qmmm_rename_residues.each{|s|
			@settings[:qmmm_rename_residues] << YAML.load(s.gsub(/^ *- */,''))
		}
		YAML.load(qmmm_rename_residues.to_s)
		@settings[:qmmm_add_ter] = '%atomname(H2)&(:BCC|:BCG)'

		# Print the output, save QM region geometry and exit unless explicitly requested
		unless @settings[:qmmm_auto_run]
			Cuby::log.puts "Automated fragmentation generated following setup:"
			Cuby::log.puts "-------------8><----------------------------------------------------------------"
			Cuby::log.puts "qmmm_core: \":#{integer_list_to_string(residues, 0, ",")}|@#{integer_list_to_string(cap_orig_atoms,1)}\""
			Cuby::log.puts "qmmm_cut_bonds:"
			Cuby::log.puts qmmm_cut_bonds
			Cuby::log.puts "qmmm_rename_residues:"
			Cuby::log.puts qmmm_rename_residues
			Cuby::log.puts "qmmm_add_ter: '%atomname(H2)&(:BCC|:BCG)'"
			Cuby::log.puts "-------------8><----------------------------------------------------------------"
			Cuby::log.puts

			@settings[:qmmm_geometry_only] = true
		end

	end
	
	def integer_list_to_string(numlist, first = 0, separator = ",", range_separator = "-") # => String
		#: Build a compact list of numbers or ranges (i.e. 1,2,5-10) from an array of
		#: integers
		return '' if numlist.size == 0
		s = ''
		array = numlist.sort
		count = nil
		array.each_index{|i|
			if i == 0
				# Start first series
				s << (array[i] + first).to_s
				count = 1
			else
				if array[i] == array[i-1] + 1
					count += 1
				else
					# Close existing series
					if count > 1
						s << range_separator
						s << (array[i-1] + first).to_s
					end
					# Start new one
					s << separator
					s << (array[i] + first).to_s
					count = 1
				end
			end
		}
		# close last series if needed
		if count > 1
			s << range_separator
			s << (array.last + first).to_s
		end

		return s
	end

end

