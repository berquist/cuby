
module ProtocolGeometry
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :tools
	PROTOCOL_DESCRIPTION = "Geometry manipulations"
	#=======================================================================
	
	require "protocols/geometry/classes/actions_atomtype.rb"
	include ActionsAtomtype
	require "protocols/geometry/classes/actions_zmat.rb"
	include ActionsZmat
	require "protocols/geometry/classes/actions_elements.rb"
	include ActionsElements
	require "protocols/geometry/classes/actions_connectors.rb"
	include ActionsConnectors
	require "protocols/geometry/classes/actions_orientation.rb"
	include ActionsOrientation

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Geometry manipulation')
	end

	def prepare_protocol
		# Read geometry
		load_geometry

		# Messages to be printed in Print
		@messages = []

		# Write or not - default is yes
		@write_geo = true

		# Perform action
		case @settings[:geometry_action]
		when :cat
			concatenate
		when :partial_pdb_update
			partial_pdb_update
		when :move_closest_contact_vdw_dist
			move_closest_contact_vdw_dist
		when :pdb_renumber
			pdb_renumber
		when :pdb_selection_as_chain
			pdb_selection_as_chain
		when :selection_to_atomlist
			selection_to_atomlist
		when :connect
			connect_connectors
		when :order_random
			order_random
		when :order_chain
			order_chain
		when :order_molecules
			order_molecules
		when :rmsd_fit
			rmsd_fit
		when :atomtype_bonds
			atomtype_bonds
		when :build_zmat
			build_zmat
		when :moments_of_inertia
			moments_of_inertia
		when :element_list
			element_list
		when :element_match
			element_match
		when :add_connector
			add_connector
		when :orient
			orient
		end

		# Save geometry
		if @write_geo
			filename = @settings[:geometry_write]
			filetype = @settings[:geometry_write_format]
			if filetype == :auto
				filetype = GeoFile.type_from_filename(filename)
				unless filetype
					Cuby::error("Failed to determine geometry file format from file name '#{filename}'")
				end
			end

			# General options
			write_options = {
				:append => @settings[:geometry_write_append],
				:extra_columns => @settings[:pdb_extra_columns]
			}

			# File-format specific options
			case filetype
			when :xyz
				write_options[:second_line] = @geometry.info[:protocol_comment] if @geometry.info[:protocol_comment]
			when :mol2
				# explicit connectivity has to be defined
				case @settings[:geometry_write_connectivity]
				when :from_radii
					conn = Connectivity.new(@geometry)
					write_options[:bondlist] = conn.bond_list
				when :from_zmat
					if @geometry.info[:z_matrix]
						write_options[:bondlist] = @geometry.info[:z_matrix].bond_list
					else
						Cuby::error("Using z-matrix connectivity was requested but the geometry does not have\na z-matrix associated.")
					end
				end
			end

			@geometry.write_file(filename, filetype, write_options)
			@messages << "Geometry written to file '#{filename}'"
		end
	end

	def run_protocol(queue)
	end

	def print_protocol
		@messages.each{|message|
			Cuby::log.puts_v(:brief,  message)
		}
		Cuby::log.puts_v(:normal, "")
	end

	def cleanup_protocol
	end

	def results
		return nil
	end

	#=======================================================================
	# Helpers
	#=======================================================================
	def must_be_pdb!(geometry, name = "")
		unless  geometry[0].properties[:pdb_atom_name]
			Cuby::error("#{name} must be a PDB.")
		end
	end

	#=======================================================================
	# Actions
	#=======================================================================
	
	def concatenate
		# Load secondary geometry
		geometry2 = Geometry.from_file(@settings[:geometry2])

		geometry2.each{|atom|
			@geometry << atom
		}
	end
	
	def partial_pdb_update
		# Load secondary geometry
		geometry2 = Geometry.from_file(@settings[:geometry2])
		# Both must be PDBs
		must_be_pdb!(@geometry, "Primary geometry")
		must_be_pdb!(geometry2, "Geometry2")

		matched_count = 0;

		@geometry.each{|atom|
			geometry2.each{|atom2|
				if atom.properties[:pdb_atom_no] == atom2.properties[:pdb_atom_no] &&
				   atom.properties[:pdb_atom_name] == atom2.properties[:pdb_atom_name] &&
				   atom.properties[:pdb_res_no] == atom2.properties[:pdb_res_no] &&
				   atom.properties[:pdb_res_name] == atom2.properties[:pdb_res_name]
					matched_count += 1;
					atom.set_coord(atom2);
				end

			}
		}

		@messages << "Coordinates of #{matched_count} atoms updated"
	end

	def move_closest_contact_vdw_dist
		# Search for two molecules
		connectivity = Connectivity.new(@geometry)
		molecules = connectivity.molecules
		Cuby::error ("Geometry / move_closest_contact_vdw_dist: Only one molecule found") if molecules.size == 1
		Cuby::error ("Geometry / move_closest_contact_vdw_dist: More than two molecules found") if molecules.size == 1
		# Get closest contact and indexes
		distance, i, j = @geometry.closest_contact(molecules[0], molecules[1])
		# Compare with vdW distance
		vdw = PeriodicTable.vdw_radius(@geometry[i].element) + PeriodicTable.vdw_radius(@geometry[j].element)
		shift = distance - vdw
		direction = (@geometry[j] - @geometry[i]) / distance * -1
		shift = direction * shift
		# Make geometry for the second molecule
		mol2 = @geometry.geometry_from_list(molecules[1])
		mol2.translate!(shift)
		@messages << "Molecules moved from distance #{'%.2f' % distance} to #{'%.2f' % @geometry[i].distance(@geometry[j])} A"
	end

	def selection_to_atomlist
		@write_geo = false # This action does not write geometry file

		atomlist = @geometry.atomlist_from_selection(@settings[:geometry_action_selection])
		@messages <<  GeometrySelections.atomlist_to_string(atomlist)
	end

	def rmsd_fit
		# Geometry2 is used as a template
		geometry2 = Geometry.from_file(@settings[:geometry2])

		@messages << "RMSD before fit: #{@geometry.rmsd_noshift(geometry2)} A"
		fitted = @geometry.fit_rmsd(geometry2)
		@messages << "RMSD after fit:  #{fitted.rmsd(geometry2)} A"
		@geometry = fitted
	end

	def moments_of_inertia
		@write_geo = false # This action does not write geometry file

		axes, eigenvalues = geometry.moment_of_inertia_diagonalized

		if @settings.set?(:geometry2)
			# When a second geometry is provided, calculate differences
			geometry2 = Geometry.from_file(@settings[:geometry2])
			axes2, eigenvalues2 = geometry2.moment_of_inertia_diagonalized

			differences = []
			eigenvalues.each_index{|i|
				differences << eigenvalues[i] - eigenvalues2[i]
			}

			@messages << "Difference of principal moments of inertia, g/mol * A^2"
			differences.each{|x| @messages << sprintf("%16.6f", x * UNIT2GMOL)}
		else
			# Just print the moments
			@messages << "Principal moments of inertia, g/mol * A^2"
			eigenvalues.each{|x| @messages << sprintf("%16.6f", x * UNIT2GMOL)}
		end
	end

	#=======================================================================
	# Actions - PDB
	#=======================================================================
	
	def pdb_renumber
		@geometry.renumber_pdb_atoms!
		@geometry.renumber_pdb_residues!
	end

	def pdb_selection_as_chain
		chain = @settings[:pdb_set_chain]
		Cuby::error "PDB chain must be one character" unless chain.length == 1
		g2 = @geometry.geometry_from_selection(@settings[:geometry_action_selection])
		g2.each{|atom|
			atom.properties[:pdb_chain] = chain
		}
		@messages << "PDB chain set to #{chain} on #{g2.size} atoms"
	end


	#=======================================================================
	# Actions - ordering
	#=======================================================================
	
	def order_random
		@geometry.shuffle!
	end

	def order_molecules
		newgeo = Geometry.new
		conn = Connectivity.new(@geometry)
		conn.molecules.each{|molecule|
			molecule.each{|i|
				newgeo << @geometry[i]
			}
		}
		@geometry = newgeo
	end

	def order_chain
		# Split geometry into terminal an backbone
		conn = Connectivity.new(@geometry)

		Cuby::error "Not implemented yet"

	end
	
end
