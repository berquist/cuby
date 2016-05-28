################################################################################
#
# Fragmentation interface
#
# Author: Jan Rezac
# Date created: 2015-09-14
# License: Cuby4 license
# Description: Fragemented calculations
# Status: Works
#
################################################################################

require "interfaces/fragmentation/classes/capped_fragments.rb"

# Missing:
# * Fragment charge definition
# * Gradient

module InterfaceFragmentation
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
		InputBlock[:calculation,          :required, "Setup for the calculation"]
	]
	#=======================================================================

	def prepare_interface
		# Build the fragments, print fragmentation info
		@fragments = CappedFragments.new(@geometry, @settings)

		# Exit if only geometry generation was requested
		if @settings[:fragmentation_geometry_only]
			if @settings[:prepare_only]
				Cuby::log.puts_debug("fragmentation_geometry_only set, skipping prepare")
				return
			else
				Cuby::log.puts_debug("fragmentation_geometry_only set, exitting")
				exit 
			end
		end

		# Set up calculations
		@calculations = []
		@fragments.all_fragments.each{|frag|
			# Debug printing
			Cuby::log.puts_debug "Preparing #{frag.id}"

			# Set up block in settings
			block_name = ("calculation" + frag.id).to_sym
			@settings.new_block(block_name)
			@settings.block(block_name).copy_from!(@settings.block(:calculation), :keep)

			# Create calculation
			calculation = Calculation.new(@name + "_" + frag.id, @settings.block(block_name), frag.geometry)
			# Link it to the geometry
			frag.calculation = calculation
			# Add to list of calculations
			@calculations << calculation
		}

		# Prepare calculations
		@calculations.each{|calc| calc.prepare(@what)}
	end

	def queue_interface(queue)
		# Update coordinates of the link atoms
		@fragments.update_link_atoms

		# Run calculations
		@calculations.each{|calculation|
			calculation.send_to_queue(queue)
		}
	end

	def compose_interface
		return composed_results
	end

	def cleanup_interface
		@calculations.each{|calculation|
			calculation.cleanup
		}
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def composed_results
		results = Results.new
		results.energy = 0.0

		# Add monomer energies
		@fragments.monomers.each{|mono|
			results.energy += mono.calculation.results.energy
		}

		# Add pairwise interaction energies
		@fragments.dimers.each_index{|i|
			i.times{|j|
				results.energy += @fragments.dimers[i][j].calculation.results.energy
				results.energy -= @fragments.monomers[i].calculation.results.energy
				results.energy -= @fragments.monomers[j].calculation.results.energy
			}
		}

		if @what.include?(:gradient)
			results.gradient = Gradient.zero(@geometry.size)
			# Monomers
			@fragments.monomers.each{|mono|
				add_gradient!(results.gradient, mono, 1.0)
			}
			# Dimers
			@fragments.dimers.each_index{|i|
				i.times{|j|
					add_gradient!(results.gradient, @fragments.dimers[i][j], 1.0)
					add_gradient!(results.gradient, @fragments.monomers[i], -1.0)
					add_gradient!(results.gradient, @fragments.monomers[j], -1.0)
				}
			}
		end
		return results
	end

	def add_gradient!(full_grad, fragment, weight)
		# Gradient from the fragment calculation
		fragment_grad = fragment.calculation.results.gradient
		# Real atoms
		fragment.atomlist_in_parent.each_index{|i|
			full_grad[fragment.atomlist_in_parent[i]] += fragment_grad[i] * weight
		}

		# Link atoms
		fragment.caps.each_index{|i|
			cap = fragment.caps[i]
			link_index = fragment.atomlist_in_parent.size + i
			g = cap.dist_ratio

			# Project gradient
			full_grad[cap.atom_from_i] += (1.0-g) * fragment_grad[link_index] * weight
			full_grad[cap.atom_orig_i] += g * fragment_grad[link_index] * weight
		}

		return nil
	end
end

