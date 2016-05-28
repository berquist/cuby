require "classes/objects/atomic_charges.rb"

class AtomicCharges
	module EEM
		# Electronegativity equalization method
		# J. Phys. Chem. A, 2002, 106 (34), pp 7887-7894
		# DOI: 10.1021/jp0205463

		# Parameters:
		# hash, indexed by element symbol, of arrays [x,n]
		# x is effective electronegativity
		# n i effective atom hardness
		# units : eV (as used in the original paper)

		# Additional parameters available in http://pubs.acs.org/doi/full/10.1021/jp020547v

		PARAMETERS_MULLIKEN = {
			# x, n
			:H	=>	[ 1.0,   17.95],
			:C	=>	[ 5.25,  9.0],
			:N	=>	[ 8.8,   9.39],
			:O	=>	[ 14.72, 14.34],	
			:F	=>	[ 15.0,  19.77]
		}

		PARAMETERS_NBO = {
			# x, n
			:H	=>	[1.0,   19.44],
			:C	=>	[8.49,  9.15],
			:N	=>	[13.45, 10.64],
			:O	=>	[27.06, 19.63],
			:F	=>	[39.18, 44.1]
		}

		def EEM.charges(geometry, charge, parameters = PARAMETERS_MULLIKEN)
			# Lin. equations matrices
			n = geometry.size
			m = Matrix.zero(n+1,n+1)
			c = Matrix.zero(n+1,1)
			n.times{|i|
				# Row
				n.times{|j|
					if i == j
						m[i,j] = 2.0 * parameters[geometry[i].element][1] * EV2KCAL * KCAL2HARTREE
					else
						m[i,j] = 1.0 / geometry[i].distance(geometry[j]) / ANGSTROM2BOHR
					end
				}
				# Close row
				m[i,n] = -1.0
				# Column matrix
				c[i,0] = -1.0 * parameters[geometry[i].element][0] * EV2KCAL * KCAL2HARTREE
			}
			# Last row: charge sum equal to system charge
			n.times{|j| m[n,j] = 1.0}
			m[n,n] = 0.0
			# Last row, column matrix
			c[n,0] = charge

			r = m.inverse * c

			charges = AtomicCharges.new(:eem)
			n.times{|i| charges << r[i,0]}
			return charges
		end
	end
end
