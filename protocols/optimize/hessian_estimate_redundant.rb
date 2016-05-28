module HessianEstimateRedundant

	def HessianEstimateRedundant.simple(internal_coords, coord_list)
		# Estimate based on coordinate types
		# Values are empirical rather than from Bakken, Helgaker, J.Chem.Phys 117(20), 2002

		v = Vector.of_size(coord_list.size)
		coord_list.each_index{|i|
			ic = internal_coords[coord_list[i]]
			case ic.type
			when :distance
				if ic.name.to_s =~ /inter/
					v[i] = 500 # 250.0
				else
					v[i] = 3200.0 # 1600.0 # 800.0
				end
			when :angle
				v[i] = 1280.0 # 640.0 # 320.0
				if ic.name == :intermolecular
					v[i] = 25.0
				end
			when :torsion
				v[i] = 640.0 # 320.0 # 160.0
				if ic.name == :intermolecular
					v[i] = 40.0 # 20.0
				end
			end
		}
		return v
	end

	def HessianEstimateRedundant.fischer_almlof(q_set, coord_list, geometry)
		#: Hessian diagonal in redundant internal coordinates, force constant
		#: derived from bond distances and table of covalent radii.
		#
		#: Reference: Fischer, Almlof, J.Phys.Chem 96(24), 9768-9774, 1992

		v = Vector.of_size(coord_list.size)

		conn = Connectivity.new(geometry)

		coord_list.each_index{|i|
			ic = q_set[coord_list[i]]
			case ic.type
			when :distance
				rcov = (PeriodicTable.covalent_radius(ic.atoms[0].element) +
					PeriodicTable.covalent_radius(ic.atoms[1].element)) *
					ANGSTROM2BOHR
				r = ic.value * ANGSTROM2BOHR
				v[i] = 0.3601 * Math.exp(-1.944 * (r - rcov))
				v[i] *=  HARTREE2KCAL / ANGSTROM2BOHR**2
			when :angle
				r1 = ic.atoms[0].distance(ic.atoms[1]) * ANGSTROM2BOHR
				r2 = ic.atoms[1].distance(ic.atoms[2]) * ANGSTROM2BOHR

				r1cov = (PeriodicTable.covalent_radius(ic.atoms[0].element) +
					PeriodicTable.covalent_radius(ic.atoms[1].element)) * ANGSTROM2BOHR
				r2cov = (PeriodicTable.covalent_radius(ic.atoms[1].element) +
					PeriodicTable.covalent_radius(ic.atoms[2].element)) * ANGSTROM2BOHR

				# It seems that there's an error in the paper in sign of ehe exponent
				# The force constants are too high unless it is corrected:
				#v[i] = 0.089 + 0.11 / (r1cov * r2cov)**-0.42 * Math::exp(-0.44*(r1 + r2 - r1cov - r2cov)) # From the paper
				v[i] = 0.089 + 0.11 / (r1cov * r2cov)**0.42 * Math::exp(-0.44*(r1 + r2 - r1cov - r2cov)) # Works better
				v[i] *=  HARTREE2KCAL
			when :torsion
				r = ic.atoms[1].distance(ic.atoms[2]) * ANGSTROM2BOHR
				rcov = (PeriodicTable.covalent_radius(ic.atoms[1].element) +
					PeriodicTable.covalent_radius(ic.atoms[2].element)) * ANGSTROM2BOHR
				l = conn.hybridization(ic.indexes[1]) + conn.hybridization(ic.indexes[2]) - 2.0

				v[i] = 0.0015 + 14.0*l**0.57/(r*rcov)**4 * Math::exp(-2.85*(r - rcov))
				v[i] *=  HARTREE2KCAL
			end
		}
		return v
	end

	def HessianEstimateRedundant.lindh(q_set, coord_list, geometry)
		#: Hessian diagonal in redundant internal coordinates, force constant
		#: derived from bond distances and constants tabulated for first 3 periods
		#: of periodic table.
		#
		#: Reference: Lindh, Bernhardsson, Karlstrom, Malmqvist, Chem.Phys.Lett 241, 423-428, 1995
		alpha_ij =    [[1.0,	0.3949,	0.3949],
			       [0.3949,	0.2800,	0.2800],
			       [0.3949,	0.2800,	0.2800]]
		rref_ij =     [[1.35,	2.1,	2.53],
			       [2.1,	2.87,	3.4],
			       [2.53,	3.4,	3.4]]
		k_dist = 0.45
		k_angle = 0.15
		k_torsion = 0.005

		v = Vector.of_size(coord_list.size)

		coord_list.each_index{|i|
			ic = q_set[coord_list[i]]
			case ic.type
			when :distance
				r = ic.value * ANGSTROM2BOHR
				p1 = PeriodicTable.period(ic.atoms[0].element) - 1
				p2 = PeriodicTable.period(ic.atoms[1].element) - 1
				p1 = 2 if p1 > 2 # Limit the indexes - JR
				p2 = 2 if p2 > 2 # Limit the indexes - JR
				v[i] = k_dist * Math.exp(alpha_ij[p1][p2]*(rref_ij[p1][p2]**2 - r**2))
				v[i] *=  HARTREE2KCAL / ANGSTROM2BOHR**2
			when :angle
				r1 = ic.atoms[0].distance(ic.atoms[1]) * ANGSTROM2BOHR
				r2 = ic.atoms[1].distance(ic.atoms[2]) * ANGSTROM2BOHR
				p1 = PeriodicTable.period(ic.atoms[0].element) - 1
				p2 = PeriodicTable.period(ic.atoms[1].element) - 1
				p3 = PeriodicTable.period(ic.atoms[2].element) - 1
				p1 = 2 if p1 > 2 # Limit the indexes - JR
				p2 = 2 if p2 > 2 # Limit the indexes - JR
				p3 = 2 if p3 > 2 # Limit the indexes - JR
				v[i] = k_angle *
					Math.exp(alpha_ij[p1][p2]*(rref_ij[p1][p2]**2 - r1**2)) *
					Math.exp(alpha_ij[p2][p3]*(rref_ij[p2][p3]**2 - r2**2))
				v[i] *=  HARTREE2KCAL
			when :torsion
				r1 = ic.atoms[0].distance(ic.atoms[1]) * ANGSTROM2BOHR
				r2 = ic.atoms[1].distance(ic.atoms[2]) * ANGSTROM2BOHR
				r3 = ic.atoms[2].distance(ic.atoms[3]) * ANGSTROM2BOHR
				p1 = PeriodicTable.period(ic.atoms[0].element) - 1
				p2 = PeriodicTable.period(ic.atoms[1].element) - 1
				p3 = PeriodicTable.period(ic.atoms[2].element) - 1
				p4 = PeriodicTable.period(ic.atoms[3].element) - 1
				p1 = 2 if p1 > 2 # Limit the indexes - JR
				p2 = 2 if p2 > 2 # Limit the indexes - JR
				p3 = 2 if p3 > 2 # Limit the indexes - JR
				p4 = 2 if p4 > 2 # Limit the indexes - JR
				v[i] = k_torsion *
					Math.exp(alpha_ij[p1][p2]*(rref_ij[p1][p2]**2 - r1**2)) *
					Math.exp(alpha_ij[p2][p3]*(rref_ij[p2][p3]**2 - r2**2)) * 
					Math.exp(alpha_ij[p3][p4]*(rref_ij[p3][p4]**2 - r3**2))
				v[i] *=  HARTREE2KCAL
			end
		}
		return v
	end

	def HessianEstimateRedundant.swart(q_set, coord_list, geometry)
		#: Hessian diagonal in redundant internal coordinates, force constant
		#: derived from bond distances
		
		# Reference: 
		
		k_dist = 0.35
		k_angle = 0.14
		k_torsion = 0.005

		v = Vector.of_size(coord_list.size)
		coord_list.each_index{|i|
			ic = q_set[coord_list[i]]
			case ic.type
			when :distance
				r1 = ic.value
				c1 = PeriodicTable.covalent_radius(ic.atoms[0].element) + PeriodicTable.covalent_radius(ic.atoms[1].element)
				s1 = Math::exp(-(r1/c1 - 1))
				v[i] = k_dist * s1
				v[i] *=  HARTREE2KCAL / ANGSTROM2BOHR**2
			when :angle
				r1 = ic.atoms[0].distance(ic.atoms[1])
				r2 = ic.atoms[1].distance(ic.atoms[2])
				c2 = PeriodicTable.covalent_radius(ic.atoms[1].element) + PeriodicTable.covalent_radius(ic.atoms[2].element)
				c1 = PeriodicTable.covalent_radius(ic.atoms[0].element) + PeriodicTable.covalent_radius(ic.atoms[1].element)
				s1 = Math::exp(-(r1/c1 - 1))
				s2 = Math::exp(-(r2/c2 - 1))
				v[i] = k_angle * s1 * s2
				v[i] *=  HARTREE2KCAL
			when :torsion
				r1 = ic.atoms[0].distance(ic.atoms[1])
				r2 = ic.atoms[1].distance(ic.atoms[2])
				r3 = ic.atoms[2].distance(ic.atoms[3])
				c1 = PeriodicTable.covalent_radius(ic.atoms[0].element) + PeriodicTable.covalent_radius(ic.atoms[1].element)
				c2 = PeriodicTable.covalent_radius(ic.atoms[1].element) + PeriodicTable.covalent_radius(ic.atoms[2].element)
				c3 = PeriodicTable.covalent_radius(ic.atoms[2].element) + PeriodicTable.covalent_radius(ic.atoms[3].element)
				s1 = Math::exp(-(r1/c1 - 1))
				s2 = Math::exp(-(r2/c2 - 1))
				s3 = Math::exp(-(r3/c3 - 1))
				v[i] = k_torsion * s1 * s2 * s3
				v[i] *=  HARTREE2KCAL
			end
		}
		return v
	end

end
