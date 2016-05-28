################################################################################
#
# Class PointCharges
#
# Author: Jan Rezac
# Date created: 2013-07-18
# License: Cuby license
# Description: Array of point charges
# Status: Works
#
################################################################################

#: Array of point charges

require "classes/objects/point_charge.rb"

class PointCharges < Array

	attr_accessor :gradient # Gradient on the point charges, if available

	#=======================================================================
	# Constructors
	#=======================================================================

	def PointCharges.from_s(string)
		#: Load a list of point charges from a string, non-empty
		#: non-comment lines are passed to PointCharge.from_s
		#: (x y z c format).

		pchs = PointCharges.new

		string.each_line{|line|
			next if line =~ /^\s*$/ # ignore blank lines
			next if line =~ /^\s*#/ # ignore comments
			pchs << PointCharge.from_s(line)
		}

		return pchs
	end

	def PointCharges.load_from_settings(string)
		#: Load point charges list from any format used in the input
		#: file. For now, only reading from string is supported.
		return PointCharges.from_s(string)
	end

	#=======================================================================
	# Printing
	#=======================================================================

	def to_g(short = false)
		#: Print the list of charges using nice formating
		return self.map{|pch| pch.to_s(short)}.join("\n")
	end

	def inspect
		return to_s
	end

	#=======================================================================
	# Calculation of interaction with molecule
	#=======================================================================
	
	def interaction_with_molecule(geometry, atomic_charges, gradient = nil, gradient_on_charges = false)
		#: Calculates energy of interaction of a molecule with given atomic charges with
		#: the curernt set of point charges. If a gradient is provided as an argument,
		#: the gradient of this interaction is added to it.
		#: When gradient_on_charges is set to true, a gradient on the point charges is
		#: also evaluated and saved within the curernt PointCharges object.

		if gradient_on_charges
			@gradient = Gradient.zero(self.size)
		end

		energy = 0.0
		geometry.each_index {|i|
			atom = geometry.at(i)
			q_atom = atomic_charges[i]
			self.each_index {|chi|
				pch = self.at(chi)
				rij = pch.distance(atom)
				energy += pch.charge * q_atom / rij / ANGSTROM2BOHR * HARTREE2KCAL
				if gradient
					df = pch.charge * q_atom * HARTREE2KCAL * -1.0 / rij**2 /  ANGSTROM2BOHR
					dcart = (pch - atom) / rij * df
					gradient[i] -= dcart
					@gradient[chi] += dcart if gradient_on_charges
				end
			}
		}

		return energy
	end

end
