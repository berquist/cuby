module ActionsAtomtype
	def atomtype_bonds
		# No geometry is written
		@write_geo = false

		# Build connectivity
		connectivity = Connectivity.new(@geometry)

		# Atom types
		atomtypes = connectivity.atom_types_simple(@geometry)

		# Traverse bonds
		connectivity.each_bond{|i,j|
			r = @geometry[i].distance(@geometry[j])
			t1, t2 = [atomtypes[i], atomtypes[j]].sort
			e1 = t1.gsub(/:.*/, '')
			e2 = t2.gsub(/:.*/, '')
			printf("%-4s%-4s%-20s%-20s%12.4f%10d%10d\n", e1,e2,t1, t2, r, i+1,j+1)
		}
	end
end
