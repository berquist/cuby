# PDB-specific methods

module GeometryPdb
	
	def is_pdb?
		return at(0).properties.has_key?(:pdb_res_name)
	end

	def residues_as_atomlists
		residues = []
		last = -1
		this_residue = nil
		each_with_index{|atom, i|
			if atom.properties[:pdb_res_no] != last
				# Create new residue
				this_residue = []
				residues << this_residue
				last = atom.properties[:pdb_res_no]
			end
			# Add the atom to the current residue
			this_residue << i
		}
		return residues
	end

	def renumber_pdb_atoms!
		 each_with_index{|atom, i|
			 atom.properties[:pdb_atom_no] = i+1
		 }
		 return nil
	end

	def renumber_pdb_residues! # => nil
		#: Fix numbering of residues in pdb. Number is increased when:
		#* original residue number changes
		#* residue name changes
		#* TER record is found
		rn = 0
		last = -1
		last_terafter = false
		last_name = ''
		each{|atom|
			if atom.properties[:pdb_res_no] != last ||
			   atom.properties[:pdb_res_name] != last_name ||
			   last_terafter
				rn += 1
			end
			last = atom.properties[:pdb_res_no]
			last_name = atom.properties[:pdb_res_name]
			last_terafter = atom.properties[:terafter]
			atom.properties[:pdb_res_no] = rn
		}
		return nil
	end

	def pdb_count_residues(strict = false) # => Integer
		#: Returns number of residues in the system
		last_res_no = first.properties[:pdb_res_no] - 1
		count = 0
		each {|atom|
			rn = atom.properties[:pdb_res_no]
			step = rn - last_res_no
			if (step < 0 || step > 1) && strict
				Cuby::error "Error in pdb file, discontinuity in numbering of residues"
			end
			if step != 0
				count += 1
			end
			last_res_no = rn
		}
		return count
	end

	def pdb_name_atoms_by_num
		#: Assign unique names to atoms. The name is created from the element abbreviation and a number (up to two digits).
		#: If there is more than 99 atoms of the same elements, error is raised. 

		numbers = atom_numbering_by_element

		# Get max. number of atoms of one element
		maxnum = numbers.max
		if maxnum > 99
			Cuby::error "Too many atoms of the elements to fit PDB numbering range 1-99"
		end

		each_index{|i|
			atom = at(i)
			e = atom.element
			atom.properties[:pdb_atom_name] = e.to_s + numbers[i].to_s
		}

		return nil
	end

end
