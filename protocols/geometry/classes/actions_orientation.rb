module ActionsOrientation

	#=======================================================================
	# Actions - orienting molecules
	#=======================================================================
	
	def orient
		# Separate selections
		selections = @settings[:geometry_orientation].strip.split(";")
		unless selections.size == 3
			Cuby::error "Keyword geometry_orientation must define three selections"
		end

		# Convert them to atomlists
		atomlists = selections.map{|s|
			list = @geometry.atomlist_from_selection(s)
			Cuby::error "Selection in geometry_orientation selects no atoms" unless list.size > 0
			list
		}

		@geometry.orient!(atomlists[0],atomlists[1],atomlists[2])
	end

end
