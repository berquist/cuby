################################################################################
#
# Interface UFF
#
# Author: Jan Rezac
# Date created: 2012-06-01
# License: Cuby4 license
# Description: Universal force filed
# Status: Partial implementation
#
################################################################################

#===============================================================================
# UFF forcefield
# J. Am. Chem. Soc. 114(25), 10024, 1992
# Partial implementation:
# * Limited set of elements
# * No electrostatics (calculation of charges difficult to implement)
#===============================================================================

require "classes/geometry/rings"
require "interfaces/uff/uff_charges"
require "classes/objects/atomic_charges_equilibration.rb"

module InterfaceUff
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :warning
	DEVELOPMENT_STATUS = "Only partial implementation of UFF, electrostatics uses different charges"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient]
	MODIFIER = false
	#=======================================================================

	# Private constants
	PARAMETER_FILE = "interfaces/uff/uff.dat"

	# Parameters class
	class UFFParameters
		attr_reader :r1, :theta0, :x1, :d1, :zeta, :z1, :vi, :uj, :xi, :hard, :radius

		def initialize(strings_array, line)
			Cuby.error("Wrong number of records on line #{line} of UFF parameter file", self) unless strings_array.size == 11
			@r1, @theta0, @x1, @d1, @zeta, @z1, @vi, @uj, @xi, @hard, @radius = strings_array.map{|x| x.to_f}
			return nil
		end
	end

	# Parameters stored in class variable
	@@uff_data = nil

	#=======================================================================
	# Interface API
	#=======================================================================

	def prepare_interface
		# Set up connectivity
		@connectivity = Connectivity.new(@geometry)
		# Assign atom types
		uff_parse_atom_types(@geometry)
		# Load forcefield parameters if needed
		uff_data_load unless @@uff_data
		# Check for missing data
		uff_check_atom_types
		# Initialize charges
		if @settings.set?(:atomic_charges_read)
			@charges = AtomicCharges.from_file(@settings[:atomic_charges_read])
		else
			#@charges = UffCharges.new(@geometry, @settings[:charge])
			# UFF charges not available, use Electronegativity Equilibration method parameterized to NBO charges
			@charges = AtomicCharges::EEM.charges(@geometry, @settings[:charge], AtomicCharges::EEM::PARAMETERS_NBO)
		end
		Cuby::log.puts_debug("Atomic charges used:")
		@geometry.each_index{|i|
			Cuby::log.puts_debug(sprintf("   %2s%10.3f",@geometry[i].element.to_s, @charges[i]))
		}
	end

	def calculate_interface
		return uff_energy_grad(@what.include?(:gradient))
	end

	#=======================================================================
	# Private methods
	#=======================================================================

	def uff_energy_grad(grad = false)
		# Initialize results
		results = Results.new
		results.gradient = Gradient.zero(@geometry.size) if grad
		results.energy = 0.0

		# Bonds
		e_bonds = 0.0
		store_r_ij = []
		@geometry.size.times{|i| store_r_ij[i] = []}
		@connectivity.each_bond{|i,j|
			at_i = @geometry.at(i)
			at_j = @geometry.at(j)
			p_i = @@uff_data[@uff_atom_types[i]]
			p_j = @@uff_data[@uff_atom_types[j]]

			r = at_i.distance(at_j)

			#!# bond order - not implemented
			bond_order = 1.0

			# Equilibrium distance
			r_ij = p_i.r1 + p_j.r1
			# bond order correction
			r_bo = -0.1332 * (r_ij) * Math::log(bond_order)
			r_ij += r_bo
			# electronegativity correction
			x_i = p_i.xi
			x_j = p_j.xi
			r_en = p_i.r1 * p_j.r1  * (x_i**0.5 - x_j**0.5)**2 /
				(x_i*p_i.r1 + x_j*p_j.r1)
			r_ij += r_en

			# Save r_ij
			store_r_ij[i][j] = r_ij

			# Morse potential depth
			d_ij = bond_order * 70.0 # kcal/mol

			# Force constant
			k_ij = 664.12 * p_i.z1 * p_j.z1 / r_ij**3

			# Energy
			alpha = (k_ij / 2.0 / d_ij)**0.5
			e = d_ij * (Math::exp(-alpha*(r-r_ij)) - 1.0)**2
			e_bonds += e

			# Gradient
			if grad
				deriv = -2.0 * alpha * d_ij * (Math::exp(-alpha*(r-r_ij)) - 1.0) * Math::exp(-alpha*(r-r_ij))
				results.gradient.add_internal_dist(@geometry, deriv, i, j)
			end
		}
		results.energy += e_bonds
		results.energy_components[:bonds] = e_bonds

		# Angles
		e_ang = 0.0
		@connectivity.each_angle{|i,j,k|
			angle = Coordinate.angle(@geometry.at(i),@geometry.at(j),@geometry.at(k))

			# Natural angle theta0
			theta0 = @@uff_data[@uff_atom_types[j]].theta0 / 180.0 * Math::PI

			# Force constant
			rij = store_r_ij[i][j] || store_r_ij[j][i] # The OR operator picks the non-nil value
			rjk = store_r_ij[j][k] || store_r_ij[k][j]

			# This is not clear in the paper, but theta0 should be used in order to make the force constant
			# a constant independent on actual angle:

			#rik = (rij**2 + rjk**2 - 2.0*rij*rjk*Math::cos(angle)) # is replaced with:
			rik = (rij**2 + rjk**2 - 2.0*rij*rjk*Math::cos(theta0))

			zi = @@uff_data[@uff_atom_types[i]].z1
			zk = @@uff_data[@uff_atom_types[k]].z1

			b = 664.12 / rij / rjk

			k_ijk = b * zi * zk / rik**5 * rij * rjk *
				(rij * rjk * (1.0 - Math::cos(theta0)**2 - rik**2 * Math::cos(theta0)))

			# Energy
			if @uff_centers[j] == :sp || @uff_centers[j] == :sp2
				# Linear, planar trigonal (supported so far) and square, octahedral
				case @uff_centers[j]
				when :sp
					n = 1.0
				when :sp2
					n = 3.0
				end
				e = k_ijk / n**2 * (1.0 - Math::cos(n*angle))

				# Gradient
				if grad
					deriv = k_ijk / n * Math::sin(n*angle)
					results.gradient.add_internal_angle(geometry, deriv, i, j, k)
				end
			else
				# Others - general formulation
				c2 = 1.0 / (4.0 * (Math::sin(theta0))**2)
				c1 = -4.0 * c2 * Math::cos(theta0)
				c0 = c2 * (2.0 * (Math::cos(theta0))**2 + 1.0)
				e = k_ijk * (c0 + c1*Math::cos(angle) + c2*Math::cos(2.0*angle))

				# Gradient
				if grad
					deriv = k_ijk * (-2.0 * Math::sin(2.0*angle)*c2 - Math::sin(angle)*c1)
					results.gradient.add_internal_angle(geometry, deriv, i, j, k)
				end
			end
			e_ang += e
		}
		results.energy += e_ang
		results.energy_components[:angles] = e_ang

		# Torsions
		e_tor = 0.0
		@connectivity.each_torsion{|i,j,k,l|
			bond_order = 1.0 #!#
			# Rules based on center hybridizations
			if @uff_centers[j] == :sp3 && @uff_centers[k] == :sp3
				n = 3.0
				f0 = 180.0
				v = (@@uff_data[@uff_atom_types[j]].vi *  @@uff_data[@uff_atom_types[k]].vi)**0.5

				if PeriodicTable.group(@geometry[j].element) == 16 && PeriodicTable.group(@geometry[k].element) == 16
					# Special case : bond between sp3 chalcogens (e.g. H2O2)
					v1 = 6.8
					v2 = 6.8
					v1 = 2.0 if @geometry[j].element == :O
					v2 = 2.0 if @geometry[k].element == :O
					v = (v1 * v2)**0.5
					n = 2.0
					f0 = 90.0
				end
			elsif @uff_centers[j] == :sp2 && @uff_centers[k] == :sp2
				n = 2.0
				f0 = 180.0
				v = 5.0 * (@@uff_data[@uff_atom_types[j]].uj *  @@uff_data[@uff_atom_types[k]].uj)**0.5 *
					(1.0 + 4.18 * Math.log(bond_order))
			elsif (@uff_centers[j] == :sp2 && @uff_centers[k] == :sp3) || (@uff_centers[j] == :sp3 && @uff_centers[k] == :sp2)
				if @uff_centers[j] == :sp2
					sp2 = j; sp3 = k
				else
					sp2 = k; sp3 = j
				end
				# General case
				n = 6.0
				f0 = 0.0
				v = 1.0
				# Special case 1 : sp3 chalcogen with single bond to sp2 atom
				if bond_order == 1.0 && PeriodicTable.group(@geometry[sp3].element) == 16
					n = 2.0
					f0 = 90.0
					v = 5.0 * (@@uff_data[@uff_atom_types[j]].uj *  @@uff_data[@uff_atom_types[k]].uj)**0.5 *
						(1.0 + 4.18 * Math.log(bond_order))
				end
				# Special case 2 : the sp2 atom is bonded to another sp2 atom
				if bond_order == 1.0
					@connectivity.bound_atoms(sp2).each{|m|
						if @uff_centers[m] == :sp2
							v = 2.0
							n = 3.0
							f0 = 180.0
							break
						end
					}
				end
			else
				# No torsion otherwise
				next
			end

			# Convert f0 to rad
			f0 = f0 * Math::PI / 180.0

			# Actual torsion angle
			f = Coordinate.torsion(@geometry.at(i),@geometry.at(j),@geometry.at(k),@geometry.at(l))

			# Energy is divided by number of torsions sharing the central bond
			torsion_n = (@connectivity.hybridization(j) - 1) * (@connectivity.hybridization(k) - 1)

			e = 0.5 * v * (1.0 - Math::cos(n*f0) * Math::cos(n*f)) / torsion_n
			e_tor += e

			# Gradient
			if grad
				deriv = 0.5 * v * n * Math::sin(f*n) * Math::cos(f0*n) / torsion_n
				results.gradient.add_internal_torsion(geometry, deriv, i, j, k, l)
			end
		}
		results.energy += e_tor
		results.energy_components[:torsions] = e_tor

		#!# Inversion

		# Nonbonded interactions
		e_el = 0.0
		e_lj = 0.0
		@connectivity.each_nonbonded_above_1_3{|i,j|
			# Lennard-Jones:
			dij = (@@uff_data[@uff_atom_types[i]].d1 *  @@uff_data[@uff_atom_types[j]].d1)**0.5
			xij = (@@uff_data[@uff_atom_types[i]].x1 *  @@uff_data[@uff_atom_types[j]].x1)**0.5

			x = @geometry.at(i).distance(@geometry.at(j))
			
			e = dij * ((xij/x)**12 - 2.0*(xij/x)**6)
			e_lj += e

			# Electrostatics: 
			e_el += @charges[i] * @charges[j] / x / ANGSTROM2BOHR * HARTREE2KCAL

			# Gradient
			if grad
				deriv = dij * ((12.0*xij**6)/x**7 - (12.0*xij**12)/x**13) # Lennard-Jones
				deriv += -1.0 * @charges[i] * @charges[j] / x**2 / ANGSTROM2BOHR * HARTREE2KCAL # Electrostatics
				results.gradient.add_internal_dist(@geometry, deriv, i, j)
			end
		}
		results.energy += e_el + e_lj
		results.energy_components[:electrostatic] = e_el
		results.energy_components[:"Lennard-Jones"] = e_lj

		return results
	end

	def uff_data_load
		filename = Cuby.install_dir + "/" + PARAMETER_FILE
		Cuby.log.puts_debug "Loading UFF data from #{filename}"
		f = File.open(filename, "r")
		@@uff_data = {}
		count = 1
		f.each{|line|
			next if line =~ /^\s*#/ # skip comments
			strs = line.strip.split(/\s+/)
			atomtype = strs.shift.to_sym
			@@uff_data[atomtype] = UFFParameters.new(strs, count)
			count += 1
		}
		f.close
	end

	def uff_check_atom_types
		@geometry.each_index{|i|
			unless @@uff_data.has_key?(@uff_atom_types[i])
				Cuby::error("UFF parameters for atom type #{@uff_atom_types[i]} (atom #{i+1}) not found", self)
			end
		}
		return nil
	end

	def uff_parse_atom_types(geo)
		Cuby.log.puts_debug "UFF atom types:"
		@uff_atom_types = []
		@uff_centers = []

		# Detect atoms in aromatic rings
		rings = Rings.new(geo)
		aromatic_atoms = []
		geo.each_index{|at_i|
			aromatic_atoms << at_i if rings.aromatic_atom?(at_i)
		}

		geo.each_with_index{|atom, i|
			# Element name, padded with "_"
			s = atom.element.to_s
			s += '_' if s.size < 2

			# Hybridization geometry (including lone pairs)
			# 1 linear (angle type linear)
			# 2 trigonal (angle type linear)
			# R resonant (angle type linear)
			# 3 tetrahedral
			# 4 planar square (angle type linear)
			# 5 trigonal bipyramidal
			# 6 octahedral (angle type linear)
			sh = '*'

			hyb = @connectivity.hybridization(i)
			batoms = @connectivity.bound_atoms(i)
			Cuby::error("UFF does not handle free atoms (atom #{i+1}, #{atom.element})", self) if hyb == 0
			if atom.element == :H
				sh = ''
				sh = 'b' if hyb == 2 # Bridge hydrogen in boranes
			elsif [:C, :N, :O].include?(atom.element)
				if hyb == 1
					sh = '*'
					if atom.element == :O
						# assume sp2 =O
						sh = '2'
					end
				else
					# Calculate average bond angle
					a = 0.0
					count = 0
					hyb.times {|k| k.times{|l|
						a += Coordinate.angle(geometry[batoms[k]], atom, geometry[batoms[l]]) * 180.0 / Math::PI
						count += 1
					}}
					a = a / count
					if a < 115.0 # sp3
						sh = '3'
						@uff_centers[i] = :sp3
					elsif a < 150.0 # sp2
						sh = '2'
						@uff_centers[i] = :sp2
					else # sp
						sh = '1'
						@uff_centers[i] = :sp
					end

					# Ring detection
					if sh == '2'
						sh = 'R' if aromatic_atoms.include?(i)
					end
				end
			elsif [:F, :Cl, :Br, :I].include?(atom.element)
				sh = ''
			else

				Cuby::error("UFF interface implements only H,C,N,O elements (#{atom.element} found)", self)
			end

			if sh == '*'
				Cuby::error("UFF interface can not determine atom type of atom #{i+1}, #{atom.element}", self)
			end

			s += sh

			@uff_atom_types[i] = s.to_sym

			Cuby.log.puts_debug "%-6d%-4s%-8s" % [i, atom.element.to_s, @uff_atom_types[i].to_s]
		}
	end
end
