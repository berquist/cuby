################################################################################
#
# Module PeriodicTable
#
# Author: Jan Rezac
# Date created: 2008-11-03
# License: Cuby license
# Description: Tables of element properties
# Extends: N/A
#
################################################################################

#: This module contains tables of some atomic properties and methods to access
#: them. 

require "classes/periodic_table/isotopes.rb"
	
module PeriodicTable
	ELEMENTS = [:X, :H, :He, :Li, :Be, :B, :C, :N, :O, :F, :Ne, :Na, :Mg,
		:Al, :Si, :P, :S, :Cl, :Ar, :K, :Ca, :Sc, :Ti, :V, :Cr, :Mn,
		:Fe, :Co, :Ni, :Cu, :Zn, :Ga, :Ge, :As, :Se, :Br, :Kr, :Rb,
		:Sr, :Y, :Zr, :Nb, :Mo, :Tc, :Ru, :Rh, :Pd, :Ag, :Cd, :In, :Sn,
		:Sb, :Te, :I, :Xe, :Cs, :Ba, :La, :Ce, :Pr, :Nd, :Pm, :Sm, :Eu,
		:Gd, :Tb, :Dy, :Ho, :Er, :Tm, :Yb, :Lu, :Hf, :Ta, :W, :Re, :Os,
		:Ir, :Pt, :Au, :Hg, :Tl, :Pb, :Bi, :Po, :At, :Rn, :Fr, :Ra,
		:Ac, :Th, :Pa, :U, :Np, :Pu, :Am, :Cm, :Bk, :Cf, :Es, :Fm, :Md,
		:No, :Lr, :Rf, :Db, :Sg, :Bh, :Hs, :Mt, :Ds, :Rg, :Uub, :Uut,
		:Uuq, :Uup, :Uuh, :Uus, :Uuo]

	#=======================================================================
	# Periodic table placement
	#=======================================================================
	
	def PeriodicTable.proton_number(element)
		#: Return proton number of an element
		pn = ELEMENTS.index(element)
		raise "Not a valid element (#{element})" if pn == nil
		return pn
	end

	def PeriodicTable.period(element)
		#: Determines period of the atom's element in the periodic table (range 1-7)
		pn = PeriodicTable::proton_number(element)
		raise "Dummy atom is not in periodic table!" if pn == 0
		return 1 if pn <= ELEMENTS.index(:He)
		return 2 if pn <= ELEMENTS.index(:Ne)
		return 3 if pn <= ELEMENTS.index(:Ar)
		return 4 if pn <= ELEMENTS.index(:Kr)
		return 5 if pn <= ELEMENTS.index(:Xe)
		return 6 if pn <= ELEMENTS.index(:Rn)
		return 7 
	end

	def PeriodicTable.group(element) # => integer or nil
		#: Return group of an element, nil for lanthanides or actinides.
		period = PeriodicTable::period(element)
		case period
		when 1
			return 1 if element == :H
			return 18 if element == :He
		when 2
			n = PeriodicTable.proton_number(element) - ELEMENTS.index(:He)
			n += 10 if n > 2
			return n
		when 3
			n = PeriodicTable.proton_number(element) - ELEMENTS.index(:Ne)
			n += 10 if n > 2
			return n
		when 4
			return PeriodicTable.proton_number(element) - ELEMENTS.index(:Ar)
		when 5
			return PeriodicTable.proton_number(element) - ELEMENTS.index(:Kr)
		when 6
			n = PeriodicTable.proton_number(element) - ELEMENTS.index(:Xe)
			return nil if (3..17).include?(n)
			n -= 14 if n > 17
			return n
		when 7
			n = PeriodicTable.proton_number(element) - ELEMENTS.index(:Rn)
			return nil if (3..17).include?(n)
			n -= 14 if n > 17
			return n
		end
		return nil
	end

	#=======================================================================
	# Data from yaml files: mass
	#=======================================================================
	
	MASS_FILE = "mass.yaml"
	@@mass_data = nil

	def PeriodicTable.mass(element)
		# Load yaml data on first use
		unless @@mass_data
			File.open(File.dirname(__FILE__) + "/" + MASS_FILE) {|f|
				@@mass_data = YAML::load(f)
			}
		end

		if m = @@mass_data[element]
			return m * GMOL2UNIT
		else
			Cuby::error("Mass not tabelated for element #{element}")
		end
	end

	#=======================================================================
	# Data from yaml files: isotopes
	#=======================================================================
	
	# Isotope mass and abundance from http://www.nist.gov/pml/data/comp.cfm
	
	ISOTOPES_FILE = "isotopes.yaml"
	@@isotopes_data = nil

	def PeriodicTable.isotopes(element)
		# Load yaml data on first use
		unless @@isotopes_data
			# Load file
			File.open(File.dirname(__FILE__) + "/" + ISOTOPES_FILE) {|f|
				@@isotopes_data = YAML::load(f)
			}
			# Convert masses to internal units
			@@isotopes_data.each_value{|isotope_list|
				isotope_list.each_value{|isotope|
					isotope.mass = isotope.mass * GMOL2UNIT if isotope.mass
				}
			}
		end

		if m = @@isotopes_data[element]
			return m
		else
			Cuby::error("Mass not tabelated for element #{element}")
		end
	end

	def PeriodicTable.common_isotope(element)
		list = PeriodicTable.isotopes(element)

		Cuby::error "No isotope data for element #{element}" if list.nil?

		case list.size
		when 0
			Cuby::error "No isotope data for element #{element}"
		when 1
			return list.values[0]
		else
			abundance = 0.0
			best = nil
			list.each_pair{|num, isotope|
				if isotope.abundance && isotope.abundance > abundance
					best = isotope
					abundance = isotope.abundance
				end
			}
			Cuby::error "Can not determine most common isotope for element #{element}, abundance data missing" if best == nil
			return best
		end
	end
	
	#=======================================================================
	# Data from yaml files: covalent radius
	#=======================================================================
	
	COVALENT_RADIUS_FILE = "covalent_radius.yaml"
	@@covalent_radius_data = nil

	def PeriodicTable.covalent_radius(element, set = :default)
		# Load yaml data on first use
		unless @@covalent_radius_data
			File.open(File.dirname(__FILE__) + "/" + COVALENT_RADIUS_FILE) {|f|
				@@covalent_radius_data = YAML::load(f)
			}
		end

		unless @@covalent_radius_data[set]
			raise "Unknown set of covalent radii requested.\nAvailable sets: #{@@covalent_radius_data.keys.join(', ')}"
		end

		if r = @@covalent_radius_data[set][element]
			return r
		else
			Cuby::error("Covalent radius not tabelated for element #{element}")
		end
	end

	#=======================================================================
	# Data from yaml files: vdw radius
	#=======================================================================
	
	VDW_RADIUS_FILE = "vdw_radius.yaml"
	@@vdw_radius_data = nil

	def PeriodicTable.vdw_radius(element, set = :default)
		# Load yaml data on first use
		unless @@vdw_radius_data
			File.open(File.dirname(__FILE__) + "/" + VDW_RADIUS_FILE) {|f|
				@@vdw_radius_data = YAML::load(f)
			}
		end

		unless @@vdw_radius_data[set]
			raise "Unknown set of vdw radii requested.\nAvailable sets: #{@@vdw_radius_data.keys.join(', ')}"
		end

		if r = @@vdw_radius_data[set][element]
			return r
		else
			Cuby::error("Covalent radius not tabelated for element #{element}")
		end
	end
	
	#=======================================================================
	# Data from yaml files: Electron configuration
	#=======================================================================

	# Name of the data file
	ELECTRON_CONFIGURATION_FILE = "electron_configuration.yaml"
	# Initialize empty data storage, it is filled on first use
	@@electron_configuration_data = nil

	# Define class used for storage of the data	
	class ElectronConfiguration
		attr_reader :reference
		attr_reader :shells
		attr_reader :orbitals
		attr_reader :occupations

		def to_s
			s = ""
			s << "[#{@reference}] " if reference
			@shells.each_index{|i|
				s << "#{@shells[i]}#{@orbitals[i]}#{@occupations[i]} "
			}
			return s.chomp
		end
	end

	def PeriodicTable.electron_configuration(element) # => ElectronConfiguration
		#: Retun electron configuration of the atom in gound state

		# Load yaml data on first use
		unless @@electron_configuration_data
			File.open(File.dirname(__FILE__) + "/" + ELECTRON_CONFIGURATION_FILE) {|f|
				@@electron_configuration_data = YAML::load(f)["electron_configuration"]
			}
		end
		# Return ElectronConfiguration
		rv = @@electron_configuration_data[element]
		Cuby::error("Electron configuration of element #{element.to_s} not found in data files") unless rv
		return rv
	end

	def PeriodicTable.valence_electrons(element, orbital_regexp = /./)
		#: Return number of valence electrons of the atom, optionally from selected orbital (s, p, d, ...)
		count = 0
		conf = PeriodicTable.electron_configuration(element)
		conf.orbitals.each_index{|i|
			count += conf.occupations[i] if conf.orbitals[i] =~ orbital_regexp
		}
		return count
	end

end
