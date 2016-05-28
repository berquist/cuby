################################################################################
#
# Module GeometrySelections
#
# Author: Jan Rezac
# Date created: 2008-11-04
# License: Cuby license
# Description: Parser of atom selection expressions
# Extends: Geometry
# Status: Tested, Documented
#
################################################################################

#: Parser of atom selection expressions. Full description of the syntax can be
#: found at the Cuby wiki.

### can be done later:
### wildchars in residue names

require "classes/misc/mod_deep_copy"
require "classes/misc/expression_parser.rb"

# Modify parser so that it works on geometry
class SelectionParser < ExpressionParser
	attr_reader :geometry
	def parse_on_geometry(geometry, string)
		@geometry = geometry
		return parse(string.gsub(/\s+/,"")).sort
	end
end

module GeometrySelections

	# Initilize the parser
	@@selection_parser = SelectionParser.new

	@@selection_parser.leaf_processor(Proc.new{|str| 
		l = @@selection_parser.geometry.atomlist_from_subselection(str)
		l = [] if l == nil
		l
	})

	@@selection_parser.add_operator("|", 1, Proc.new{|l,r| l | r})
	@@selection_parser.add_operator("&", 0, Proc.new{|l,r| l & r})

	@@selection_parser.add_function("%atomname", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_atomname(parameters)
		})

	@@selection_parser.add_function("%pdb_no", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_pdb_no(parameters)
		})

	@@selection_parser.add_function("%coord", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_coord(parameters)
		})

	@@selection_parser.add_function("%same_residue", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_same_residue(parameters)
		})

	@@selection_parser.add_function("%within", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_within(parameters)
		})

	@@selection_parser.add_function("%all", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_all(parameters)
		})

	@@selection_parser.add_function("%molecule", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_molecule(parameters)
		})

	@@selection_parser.add_function("%not", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_not(parameters)
		})

	@@selection_parser.add_function("%pdb_chain", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_pdb_chain(parameters)
		})

	@@selection_parser.add_function("%pdb_segid", 
		Proc.new{|parameters| 
			@@selection_parser.geometry.selection_func_pdb_segid(parameters)
		})



	#=======================================================================
	# The almighty method and one derivative of it
	#=======================================================================
	
	def atomlist_from_selection(selection)
		#: Get list of atoms (array of indexes) that match the selection

		# Reset connectivity
		@connectivity = nil
		@molecules = nil

		# Parse expression using RDParser
		atomlist = @@selection_parser.parse_on_geometry(self, selection)
	end

	def geometry_from_selection(selection)
		#: Get Geometry object of atoms that match the selection
		return geometry_from_list(atomlist_from_selection(selection)) 
	end

	def geometry_from_list(atomlist)
		g = Geometry.new
		atomlist.each{|i|
			g << at(i)
		}
		return g
	end

	#=======================================================================
	# Set operations
	#=======================================================================
	
	def GeometrySelections.overlap?(atomlist1, atomlist2) # => true/false
		#: Check overlap of two atomlists
		return (atomlist1 & atomlist2).size > 0
	end

	#=======================================================================
	# Reverse: string from atomlist
	#=======================================================================

	def GeometrySelections.atomlist_to_string(atomlist, first = 1) # => String
		#: Convert atom list (array of indexes) to selection expression
		#: in the simple syntax (list of index ranges). 
		separator = ","
		range_separator = "-"
		return '' if atomlist.size == 0
		s = ''
		array = atomlist.sort
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

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def atomlist_from_subselection(selection)
		#: Parse elementary selection expression
		
		# check for allowed characters
		if selection =~ /[^a-z,A-Z,0-9,~,\-,\,,:,@,*]/
			Cuby::error("Selection \"#{selection}\" contains illegal characters")
		end

		# check if it is simple selection ([0-9,\,,-] only)
		if selection =~ /^[0-9,\,,-]+$/
			return atomlist_from_simple_selection(selection)
		end

		# look for components
		has_residues = selection =~ /:/
		has_atoms    = selection =~ /@/

		residue_part = selection.gsub(/@.*/,"").gsub(/^:/,"")
		atom_part    = selection.gsub(/.*@/,"")

		if has_residues
			# check if residues exist in geometry
			if at(0).properties[:pdb_res_name] == nil
				Cuby::warning("Selection of residues can be done only on data from PDB files, nothing will be selected", :once)
				each {|atom|
					atom.properties[:pdb_res_name] = 'UNK'
					atom.properties[:pdb_res_no] = 1
				}
			end

			# handle negation operator
			negate_res = false
			if residue_part =~ /^~/
				negate_res = true
				residue_part.gsub!(/^~/,'')
			end

			# check for acceptable characters again (no wildcards here)
			Cuby::error("Selection of residues in \"#{selection}\" contains illegal characters") unless residue_part =~ /^[A-Z,a-z,0-9,\,,\-]+$/
			# split by ","
			res_sel_pieces = residue_part.split(',')
			reslist = res_sel_pieces.collect{|piece| residuelist_from_residue_selection_fragment(piece)}.flatten.uniq

			# handle negation
			if negate_res
				# build list of all residues (by number)
				allres = []
				each {|atom|
					allres << atom.properties[:pdb_res_no] unless allres.include?(atom.properties[:pdb_res_no])
				}
				reslist = allres - reslist
			end

			# create atomlist from reslist
			atomlist_res = atomlist_from_reslist(reslist)
		end

		if has_atoms
			# handle negation operator
			negate_atoms = false
			if atom_part =~ /^~/
				negate_atoms = true
				atom_part.gsub!(/^~/,'')
			end

			# check for acceptable characters again (no wildcards here)
			Cuby::error("Selection of atoms in \"#{selection}\" contains illegal characters") unless atom_part =~ /^[A-Z,a-z,0-9,\,,\-]+$/

			# split by ","
			atom_sel_pieces = atom_part.split(',')
			atomlist_at = atom_sel_pieces.collect{|piece| atomlist_from_atom_selection_fragment(piece)}.flatten.uniq

			# handle negation
			if negate_atoms
				atomlist_at = (0..(self.size - 1)).to_a - atomlist_at
			end
		end

		# return the right list
		return atomlist_res if has_residues && ! has_atoms
		return atomlist_at  if ! has_residues && has_atoms
		return []  if ! has_residues && ! has_atoms
		# create final list as overlap between atomlist_at and atomlist_res
		return atomlist_res & atomlist_at # union of the lists

	end

	def atomlist_from_simple_selection(selection)
		#: Parse elementary expression in the simple syntax
		atom_sel_pieces = selection.split(',')
		atomlist_at = atom_sel_pieces.collect{|piece| atomlist_from_atom_selection_fragment(piece)}.flatten.uniq
		return atomlist_at
	end

	def atomlist_from_atom_selection_fragment(fragment)
		#: return array of indexes of atoms that matches following selections:
		#* atom number (starting from 1) (i.e. 5)
		#* atom number range (i.e. 1-10)
		#* element (all atoms of selected element)
		if fragment =~ /^[a-z,A-Z]+$/
			# it's an element
			begin
				element = Atom.element_from_string(fragment)
			rescue
				Cuby::error("Invalid fragment of atom selection \"#{fragment}\", unknown element")
			end
			return atomlist_of_element(element)
		elsif fragment =~ /^[0-9]+$/
			# it's a number
			number = fragment.to_i
			Cuby::error("Invalid fragment of atom selection \"#{fragment}\", atom numbering starts at 1") if number <= 0
			Cuby::error("Invalid fragment of atom selection \"#{fragment}\", selection beyond size of the system") if number > self.size
			return [number-1]
		elsif fragment =~ /^[0-9]+-[0-9]+$/
			# it's a range
			bounds = fragment.split('-').map{|x| x.to_i}
			start = bounds.min
			stop = bounds.max
			Cuby::error("Invalid fragment of atom selection \"#{fragment}\", atom numbering starts at 1") if start <= 0
			Cuby::error("Invalid fragment of atom selection \"#{fragment}\", selection beyond size of the system") if start > self.size
			Cuby::error("Invalid fragment of atom selection \"#{fragment}\", atom numbering starts at 1") if stop <= 0
			Cuby::error("Invalid fragment of atom selection \"#{fragment}\", selection beyond size of the system") if stop > self.size
			return ((start-1)..(stop-1)).to_a
		else
			Cuby::error("Fragment of atom selection \"#{fragment}\" should be element, number or range of numbers (i.e. 1-5)")
		end
	end

	def atomlist_of_element(element)
		#: return array of indexes of atoms of the specified element
		list = []
		each_index{|i|
			list << i if at(i).element == element
		}
		return list
	end

	def residuelist_from_residue_selection_fragment(fragment)
		#: return array of indexes of residues that matches following selections:
		#* residue number (starting from 1) (i.e. 5)
		#* residue number range (i.e. 1-10)
		#* residue name (all residues of this name)
		if fragment =~ /[a-z,A-Z][a-z,A-Z,0-9]*$/
			# it's a name
			list = residuelist_of_name(fragment)
			#Cuby::warning("Fragment of residue selection \"#{fragment}\" does not match any residues in the geometry.") if list.size == 0
			return list
		elsif fragment =~ /^[0-9]+$/
			# it's a number
			number = fragment.to_i
			Cuby::error("Invalid fragment of residue selection \"#{fragment}\", residue numbering starts at 1") if number <= 0
			return [number]
		elsif fragment =~ /^[0-9]+-[0-9]+$/
			# it's a range
			bounds = fragment.split('-').map{|x| x.to_i}
			start = bounds.min
			stop = bounds.max
			Cuby::error("Invalid fragment of residue selection \"#{fragment}\", residue numbering starts at 1") if start <= 0
			Cuby::error("Invalid fragment of residue selection \"#{fragment}\", residue numbering starts at 1") if stop <= 0
			return ((start)..(stop)).to_a
		else
			Cuby::error("Fragment of residue selection \"#{fragment}\" should be residue name, number or range of numbers (i.e. 1-5)")
		end
	end

	def residuelist_of_name(name)
		#: Get list of residue numbers for residues of given name
		reslist = []
		each{|atom|
			if atom.properties[:pdb_res_name] == name
				reslist << atom.properties[:pdb_res_no] unless reslist.include?(atom.properties[:pdb_res_no])
			end
		}
		# Check for wrong numbering
		each {|atom|
			if reslist.include?(atom.properties[:pdb_res_no])
				if atom.properties[:pdb_res_name] != name
					Cuby::error("Wrong numbering of residues in PDB,\nresidue #{name} has the same number (#{atom.properties[:pdb_res_no]}) as residue #{atom.properties[:pdb_res_name]}")
				end
			end
		}
		return reslist
	end

	def atomlist_from_reslist(reslist)
		#: Get atom list from list of residue numbers
		atomlist = []
		each_index{|i|
			atomlist << i if reslist.include?(at(i).properties[:pdb_res_no])
		}
		return atomlist
	end

	#=======================================================================
	# Functions
	#=======================================================================
	
	def selection_func_atomname(names)
		list = []
		each_index{|i|
			names.split(",").each{|name|
				list << i if at(i).properties[:pdb_atom_name] == name
			}
		}
		return list
	end

	def selection_func_pdb_no(names)
		if at(0).properties[:pdb_atom_no] == nil
			Cuby::warning("Selection of PDB atom numbers can be done only on data from PDB files, nothing will be selected", :once)
		end

		list = []
		names.split(",").each{|range|
			if range =~ /^[0-9]+$/
				start = range.to_i
				stop = range.to_i
			elsif range =~ /^[0-9]+-[0-9]+$/
				bounds = fragment.split('-').map{|x| x.to_i}
				start = bounds.min
				stop = bounds.max
			else
				Cuby::error("Invalid range specifiaction in %pdb_no() selection");
			end
			each_index{|i|
				list << i if at(i).properties[:pdb_atom_no] >= start && at(i).properties[:pdb_atom_no] <= stop
			}
		}
		list.uniq!
		return list
	end


	def selection_func_coord(inner)
		if matchdata = /([x,y,z])([>,<,=])(-?[0-9]*\.?[0-9]+)/.match(inner)
			coordinate = ["x","y","z"].index(matchdata[1])
			operator = matchdata[2]
			value = matchdata[3].to_f
			list = []
			case operator
			when "="
				each_index{|i|
					list << i if at(i)[coordinate] == value
				}
			when ">"
				each_index{|i|
					list << i if at(i)[coordinate] > value
				}
			when "<"
				each_index{|i|
					list << i if at(i)[coordinate] < value
				}
			end
		else
			Cuby::error "Invalid selection expression '#{selection}'"
		end
		return list
	end

	def selection_func_same_residue(inner)
		atomlist = atomlist_from_selection(inner)
		reslist = []
		atomlist.each{|i|
			reslist << at(i).properties[:pdb_res_no]
		}
		reslist.uniq!
		return atomlist_from_selection(":"+reslist.join(","))
	end

	def selection_func_within(parms)
		distance, inner = parms.split(';') ### here, only first semicolon shoul split the string
		sellist = atomlist_from_selection(inner)
		distance = distance.to_f
		list = []
		sellist.each{|i|
			a = at(i)
			each_index{|j|
				if a.distance(at(j)) <= distance
					list << j
				end
			}
		}
		list.uniq!
		return list
	end

	def selection_func_all(parms)
		list = (0..(size-1)).to_a
		return list
	end

	def selection_func_molecule(parms)
		#Evaluate conenctivity and list of molecules
		unless @connectivity
			@connectivity = Connectivity.new(self)
		end
		unless @molecules
			@molecules = @connectivity.molecules
		end

		mol_num, req_mols = parms.split(';')
		mol_num = mol_num.to_i

		if mol_num > @molecules.size
			Cuby::error "Only #{@molecules.size} molecules found, can't select molecule no. #{mol_num}"
		end

		if req_mols
			unless @molecules.size == req_mols.to_i
				Cuby::error "#{req_mols.to_i} molecules expected but #{@molecules.size} found"
			end
		end

		return @molecules[mol_num - 1]
	end

	def selection_func_not(parms)
		list = (0..(size-1)).to_a
		sel = atomlist_from_selection(parms)
		return list - sel
	end

	def selection_func_pdb_chain(name)
		if at(0).properties[:pdb_atom_no] == nil
			Cuby::warning("Selection of PDB chain can be done only on data from PDB files, nothing will be selected", :once)
		end

		list = []
		each_index{|i|
			list << i if at(i).properties[:pdb_chain] == name
		}
		list.uniq!
		return list
	end

	def selection_func_pdb_segid(name)
		if at(0).properties[:pdb_atom_no] == nil
			Cuby::warning("Selection of PDB segment id can be done only on data from PDB files, nothing will be selected", :once)
		end

		list = []
		each_index{|i|
			next unless at(i).properties[:pdb_segid]
			list << i if at(i).properties[:pdb_segid].strip == name
		}
		list.uniq!
		return list
	end

end
