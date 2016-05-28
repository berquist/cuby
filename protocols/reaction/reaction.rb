
module ProtocolReaction
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :composite_calculation
	PROTOCOL_DESCRIPTION = "Reaction energy calculations"
	# Input structure
	INPUT_BLOCKS = [
		InputBlock[:calculation, :required, "Setup for the calculation"]
	]
	#=======================================================================
	
	require "protocols/reaction/classes/chemical_reaction.rb"

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Reaction energy calculation')
	end

	def prepare_protocol
		# Parse reaction
		@reaction = ChemicalReaction.from_string(@settings[:reaction_formula])

		if @settings[:reaction_print].include?(:formula)
			Cuby.log.puts_v(:normal, "Reaction:")
			Cuby.log.puts_v(:normal, @reaction.to_s)
			Cuby.log.puts_v(:normal, "")
		end

		# Create child jobs
		@jobs = []
		@reaction.all_items.each{|item|
			# Create settings block
			blockname = ("calculation_" + item.name.downcase).to_sym
			@settings.new_block(blockname) unless @settings.has_block?(blockname)
			@settings.block(blockname).copy_from!(@settings.block(:calculation), :keep)
			# Create settings
			settings = @settings.block(blockname)

			# Get geometry
			if @settings[:reaction_smiles]
				settings[:geometry] = "smiles:" + item.name
			else
				settings[:geometry] = @settings[:reaction_geometries][item.name]
			end

			# Create and save the job
			newjob =Job.new(@name + "_" + item.name, settings)
			item.job = newjob
			@jobs << newjob
		}

	
		@jobs.each{|job| job.prepare}
	end

	def run_protocol(queue)
		# Submit jobs to queue
		@jobs.each{|job| job.run}
	end

	def cleanup_protocol
		# Clean child jobs
		@jobs.each{|job| job.cleanup}
	end

	def print_protocol
		results.print_energy("Reaction energy")
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			@results = Results.new
			@results.energy = 0.0
			@reaction.reactants.each{|item|
				@results.energy -= item.job.results.energy * item.num
			}
			@reaction.products.each{|item|
				@results.energy += item.job.results.energy * item.num
			}
		end
		return @results
	end

	#=======================================================================
	# Optional protocol functionality
	#=======================================================================
	
	def single_line_result
		return sprintf("Reaction energy: %20.6f kcal/mol", results.energy)
	end

	def single_number_result
		return results.energy
	end
end
