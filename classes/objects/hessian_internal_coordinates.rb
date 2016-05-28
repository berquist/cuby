
class HessianInternalCoordinates
	# Hessian in (redundant) internal coordinates. For now, only the diagonal elements are considered.
	
	attr_accessor :diag_elem_array
	
	#=======================================================================
	# DiagonalElement class
	#=======================================================================

	class DiagonalElement
		attr_accessor :atoms_types
		attr_accessor :coord_type
		attr_accessor :coord_value
		attr_accessor :hess_value
		attr_accessor :weight   # number of occurences of given atoms_types
		
		def initialize(coord_type, coord_value, atoms_types, hess_value, weight = 1)
			@coord_type = coord_type
			@coord_value = coord_value
			@atoms_types = atoms_types
			@hess_value = hess_value
			@weight = weight
		end
		
		def DiagonalElement.from_molecule(geometry, transformation, all_atoms_types, int_hessian, i)
			atoms_types = []
			transformation.int_coords[i].indexes.each do |atom_index|
				atoms_types << all_atoms_types[atom_index]
			end
			if atoms_types.first > atoms_types.last
				atoms_types.reverse!
			elsif atoms_types.size >= 4 and atoms_types.first == atoms_types.last and atoms_types[2] > atoms_types[3]
				atoms_types.reverse!
			end
			
			if int_hessian
				dement = DiagonalElement.new(transformation.int_coords[i].type, transformation.int_coords[i].value, atoms_types, int_hessian[i,i])
			else
				dement = DiagonalElement.new(transformation.int_coords[i].type, transformation.int_coords[i].value, atoms_types, 0.0)
			end
			return dement
		end
		
		def atoms_self   # returns coordinate's atoms only, without their types/neighbors
			a = []
			@atoms_types.each do |atomtype|
				a << atomtype.split(":")[0]   # returns prefix until ":", i.e. "Na:Cl" => "Na"
			end
			return a
		end
	end
	
	#=======================================================================
	# HessianInternalCoordinates methods
	#=======================================================================
	
	def initialize
		@diag_elem_array = []   # array of all diagonal elements of class DiagonalElement
	end
	
	def HessianInternalCoordinates.from_molecule(geometry, transformation, int_hessian)
		hic = HessianInternalCoordinates.new
		
		conn = Connectivity.new(geometry)
		all_atoms_types = conn.atom_types_simple(geometry)
		
		int_hessian.m.times do |i|
			hic.diag_elem_array << DiagonalElement.from_molecule(geometry, transformation, all_atoms_types, int_hessian, i)
		end
		
		return hic
	end
	
	def HessianInternalCoordinates.from_yaml(infile)
		f = File.open(infile, "r")
		hic = HessianInternalCoordinates.new
		hic.diag_elem_array = YAML.load(f)
		f.close
		return hic
	end
	
	def best_match_value(diag_element, stat_counters = nil)   # usage: create empty DiagonalElement, set only coord_type and atoms_types
		atoms_self_vals = []
		coord_type_vals = []
		
		@diag_elem_array.each do |db_item|
			if diag_element.atoms_types == db_item.atoms_types    # in case of perfect accordance (atoms + their types)
				stat_counters[:diag_hessian_perfect] += 1 if stat_counters
				return db_item.hess_value
			end
			if diag_element.atoms_self == db_item.atoms_self   # prepare a field of values for the case of imperfect accordance
				atoms_self_vals << db_item.hess_value
			end
			if diag_element.coord_type == db_item.coord_type   # prepare a field of values for the case of no accordance
				coord_type_vals << db_item.hess_value
			end
		end
		unless atoms_self_vals.empty? # if there is nonempty field of hessian values for given atoms themself, returns average
			Cuby::log.puts_debug("almost      #{diag_element.coord_type}   #{diag_element.atoms_types}")
			stat_counters[:diag_hessian_almost] += 1 if stat_counters
			return atoms_self_vals.reduce(:+).to_f / atoms_self_vals.size
		end
		unless coord_type_vals.empty?   # if there is no coordinate with these atoms
			Cuby::log.puts_v(:debug, "not found   #{diag_element.coord_type}   #{diag_element.atoms_types}")
			stat_counters[:diag_hessian_not_found] += 1 if stat_counters
			return coord_type_vals.reduce(:+).to_f / coord_type_vals.size
		end
		return 100   # hypothetically ...
	end
	
	def write_yaml(outfile)
		File.open(outfile, "w+"){|f| f.puts @diag_elem_array.to_yaml}
	end
	
	def include_by_type(dement)
		@diag_elem_array.each_index do |i|
			return [true, i] if dement.atoms_types == @diag_elem_array[i].atoms_types
		end
		return [false, -1]
	end
	
end
