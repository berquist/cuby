################################################################################
#
# Restraints interface
#
# Author: Jan Rezac
# Date created: 2012-04-11
# License: Cuby4 license
# Description: Harmonic restraints in internal coordinates
# Status: Works
#
################################################################################

require "classes/internal_coordinates/internal_coordinate.rb"

module InterfaceRestraints
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "All cuby3 functionality implemented, examples for all features"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient]
	IGNORE = [:point_charges, :point_charges_gradient]
	MODIFIER = true
	#=======================================================================


	def prepare_interface
		#: Parse settings, build lists of coordinates
		@restraints = []
		@settings[:restraints].each { |line|
			@restraints << Restraint.from_string(line, @geometry)
		}
	end

	def calculate_interface
		results = Results.new

		#: Calculate
		if what.include?(:energy)
			results.energy = restraints_energy
		end
		if what.include?(:gradient)
			results.gradient = restraints_gradient
		end

		return results
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def restraints_energy
		#: Returns sum of energies of all restraints of the system
		energy = 0.0
		@restraints.each{|res|
			energy += res.energy
		}
		return energy
	end

	def restraints_gradient
		#: Returns cartesian gradients of all restraints of the system
		gradient = Gradient.zero(@geometry.size)

		# Iterate over restraints, add the gradients to the array
		@restraints.each{|res|
			res.restraint_grad_to_grad_array(gradient)
		}
		return gradient
	end


	#=======================================================================
	# Restraint class
	#=======================================================================
	
	class Restraint < GroupsInternalCoord
		#: Class based on InternalCoord, added is reference value and force
		#: constant that is applied to reach it. 

		def Restraint.from_string(string, geometry) # => Restraint
			#: Creates new instance of Restraint from string of format used in
			#: cuby input.

			string = string.gsub(/ +/,'')
			entries = string.split(';')
			case entries[0].downcase
			when 'distance'
				Cuby::error "Wrong number of parameters for distance restraint" if entries.size != 5
				return Restraint.new(:distance, entries[1].to_f, entries[2].to_f, 0.0, geometry,
						     geometry.atomlist_from_selection(entries[3]),
						     geometry.atomlist_from_selection(entries[4]))
			when 'distance_f'
				Cuby::error "Wrong number of parameters for flat-bottom distance restraint" if entries.size != 6
				return Restraint.new(:distance, entries[1].to_f, entries[2].to_f, entries[3].to_f,  geometry,
						     geometry.atomlist_from_selection(entries[4]),
						     geometry.atomlist_from_selection(entries[5]))
			when 'angle'
				Cuby::error "Wrong number of parameters for angle restraint" if entries.size != 6
				return Restraint.new(:angle, entries[1].to_f / 180.0 * Math::PI, entries[2].to_f, 0.0, geometry,
						     geometry.atomlist_from_selection(entries[3]),
						     geometry.atomlist_from_selection(entries[4]),
						     geometry.atomlist_from_selection(entries[5]))
			when 'torsion'
				Cuby::error "Wrong number of parameters for torsion restraint" if entries.size != 7
				return Restraint.new(:torsion, entries[1].to_f / 180.0 * Math::PI, entries[2].to_f, 0.0, geometry,
						     geometry.atomlist_from_selection(entries[3]),
						     geometry.atomlist_from_selection(entries[4]),
						     geometry.atomlist_from_selection(entries[5]),
						     geometry.atomlist_from_selection(entries[6]))
			when 'distance_difference'
				Cuby::error "Wrong number of parameters for distance difference restraint" if entries.size != 7
				return Restraint.new(:distance_difference, entries[1].to_f, entries[2].to_f, 0.0, geometry,
						     geometry.atomlist_from_selection(entries[3]),
						     geometry.atomlist_from_selection(entries[4]),
						     geometry.atomlist_from_selection(entries[5]),
						     geometry.atomlist_from_selection(entries[6]))
			else
				Cuby::error "Unknown type of internal coordinate to be restrained: \"#{entries[0]}\""
			end
		end

		def initialize(type, ref_value, fconst, x0, geometry, *lists) # => Restraint
			#: Extends InternalCoord's initialize with new parameters
			super(type, geometry, *lists)
			@fconst = fconst
			@ref_value = ref_value
			@x0 = x0
		end

		def energy # => Float
			#: Returns energy of the restraint in current geometry
			dx = value - @ref_value
			if @type == :torsion
				dx += 2.0 * Math::PI if dx <= -Math::PI
				dx -= 2.0 * Math::PI if dx > Math::PI
			end
			if dx.abs <= @x0
				return 0.0
			else
				return 0.5 * @fconst * (dx-@x0)**2
			end
		end

		def internal_gradient # => Float
			#: Returns gradient of the restrain in the respective internal coordinate
			dx = value - @ref_value
			if @type == :torsion
				dx += 2.0 * Math::PI if dx <= -Math::PI
				dx -= 2.0 * Math::PI if dx > Math::PI
			end
			if dx.abs <= @x0
				return 0.0
			else
				if dx > 0
					return  @fconst * (dx-@x0)
				else
					return  @fconst * (dx+@x0)
				end
			end
		end

		def restraint_grad_to_grad_array(gradarray) # => nil
			#: Adds gradient of the restraint to array containing cartesian gradients in the geometry 
			grad_to_grad_array(internal_gradient, gradarray)
			return nil
		end

	end
end
