################################################################################
#
# Module XBond
#
# Author: J. Rezac
# Date created: 2014-06-13
# License: Cuby license
# Description: Halogen bond correction for PM6
# Extends: Calculation
# Status: Not documented, in progress
#
################################################################################

require "classes/math/func.rb"

module InterfaceXBond
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "Finished and tested"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient]
	MODIFIER = true
	#=======================================================================
	
	#=======================================================================
	# Data
	#=======================================================================
	DEFAULT_PARAMETERS = {
		:Cl => {
			:N => [1048900000000.0, -9.946],
			:O => [1871000000.00, -7.44]
		},
		:Br => {
			:N => [55600.0, -3.04],
			:O => [21600.0, -3.30],
	#		:C => [10000.0, -3.0]
		},
		:I  => {
			:N => [523744000.0, -6.77],
			:O => [2436000.0, -4.71],
			:S => [1051000.0, -3.82]
		}
	}

	#=======================================================================
	# Interface methods
	#=======================================================================

	def prepare_interface
		# Check if parent is PM6
		if @settings.parent && @settings.parent.set?(:method)
			method = @settings.parent[:method]
			if method != :pm6
				#Cuby::error "The X-bond correction is intended for use with PM6 method only,\nmethod #{method.to_s.upcase} used in the input"
			end
		end

		x_bond_load_parameters

		x_bond_make_lists
	end

	def calculate_interface
		# Reload parameters during parametrization
		if @settings[:development][:x_bond_mod_pair]
			x_bond_load_parameters
		end

		# Calculate
		results = x_bond_energy_grad(@what.include?(:gradient))
		return results
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def x_bond_load_parameters
		# Load default parameters
		@parameters = DEFAULT_PARAMETERS

		# Add parameters from input
		if @settings.set?(:xbond_parameters)
			@settings.elements_hash(:xbond_parameters).each_pair{|halogen,acceptors|
				Cuby::error "Keyword xbond_parameters should be hash of hashes" unless acceptors.class == Hash
				acceptors.each_pair{|acceptor, values|
					begin
						acceptor = Atom.element_from_string(acceptor)
					rescue
						Cuby::error "The xbond_parameters keyword should be indexed by elements only\nfound index #{acceptor}"
					end
					if values.class == Array # && values.size == 2 
						Cuby::log.puts_debug "X_bond: setting parameters for #{halogen}-#{acceptor}"
						@parameters[halogen] = {} unless @parameters[halogen]
						@parameters[halogen][acceptor] = values
					else
						Cuby::error "The values of xbond_parameters keyword should be arrays of size 2"
					end
				}
			}
		end

		if @settings[:development][:x_bond_mod_pair]
			ele1, ele2 =  @settings[:development][:x_bond_mod_pair].strip.split
			ele1 = Atom.element_from_string(ele1)
			ele2 = Atom.element_from_string(ele2)
			Cuby::log.puts_debug "Development x_bond parameters for #{ele1} - #{ele2}"
			@parameters[ele1] = {} unless @parameters[ele1]
			@parameters[ele1][ele2] = [@settings[:x0], @settings[:x1]]
			@parameters[ele1][ele2] << @settings[:x2] if @settings.set?(:x2)
			@parameters[ele1][ele2] << @settings[:x3] if @settings.set?(:x3)
			@parameters[ele1][ele2] << @settings[:x4] if @settings.set?(:x4)
			@parameters[ele1][ele2] << @settings[:x5] if @settings.set?(:x5)
			@parameters[ele1][ele2] << @settings[:x6] if @settings.set?(:x6)
		end
	end

	def x_bond_make_lists
		@pairs = [] # list of arrays [halogen, acceptor]

		if @settings[:xbond_only_noncovalent]
			connectivity = Connectivity.new(@geometry)
			connectivity.each_nonbonded_above_1_4{|i,j|
				atom_i = @geometry.at(i)
				atom_j = @geometry.at(j)
				if @parameters.has_key?(atom_i.element) && @parameters[atom_i.element].has_key?(atom_j.element)
					@pairs << [i,j]
				elsif @parameters.has_key?(atom_j.element) && @parameters[atom_j.element].has_key?(atom_i.element)
					@pairs << [j,i]
				end
			}
		else
			@geometry.size.times{|i| i.times{|j|
				atom_i = @geometry.at(i)
				atom_j = @geometry.at(j)
				if @parameters.has_key?(atom_i.element) && @parameters[atom_i.element].has_key?(atom_j.element)
					@pairs << [i,j]
				elsif @parameters.has_key?(atom_j.element) && @parameters[atom_j.element].has_key?(atom_i.element)
					@pairs << [j,i]
				end
			}}
		end
	end
	

	def x_bond_energy_grad(grad = false)
		# Initialize energy and gradient wit zeroes
		results = Results.new
		energy = 0.0
		results.gradient = Gradient.zero(@geometry.size) if grad

		@pairs.each{|array|
			halogen, acceptor = array
			at_h = @geometry.at(halogen)
			at_a = @geometry.at(acceptor)

			if para = @parameters[at_h.element][at_a.element]
				r = at_h.distance(at_a)
				# Smooth transition to flat potential occurs between 0.7 and 0.8 * rVDW
				rvdw = PeriodicTable.vdw_radius(at_h.element) + PeriodicTable.vdw_radius(at_a.element)
				r0 = rvdw * 0.7
				r1 = rvdw * 0.8
				swf = Func::PolySwitching.new(r0, r1, :sw_0_to_1)
				swf_r = swf.calc(r)

				if para.size == 4
					k, e, p, rx = para
					if r > rx
						exp = k * Math::exp(e * (r - rx)**p)
					else
						exp = k
					end
					if r0 > rx
						exp0 = k * Math::exp(e * (r0 - rx)**p)
					else
						exp0 = k
					end
				elsif para.size == 3
					k, e, p = para
					exp = k * Math::exp(e * r**p)
					exp0 = k * Math::exp(e * r0**p)
				elsif para.size == 5
					k, e, c, a,x = para
					exp = k * Math::exp(e * r)
					exp0 = k * Math::exp(e * r0)

					# Attraction
					exp -= (c/r)**a
					exp0 -= (c/r0)**a
				elsif para.size == 6
					k, e, p, c, a,x = para
					exp = k * Math::exp(e * r**p)
					exp0 = k * Math::exp(e * r0**p)

					# Attraction
					exp -= (c/r)**a
					exp0 -= (c/r0)**a
				else
					k, e = para
					exp = k * Math::exp(e * r)
					exp0 = k * Math::exp(e * r0)
				end

				energy += exp * swf_r + (1.0-swf_r) * exp0

				# Gradient
				if grad
					if para.size == 4
						deriv_e = e * k * p * (r-rx)**(p-1) * Math::exp(e * (r-rx)**p)
					elsif para.size == 3
						deriv_e = e * k * p * r**(p-1) * Math::exp(e*r**p)
					else
						deriv_e =  e * k * Math::exp(e*r)
					end
					deriv_s = swf.deriv(r)
					deriv = deriv_e * swf_r + deriv_s * exp +
						-deriv_s * exp0
					results.gradient.add_internal_dist(@geometry, deriv, halogen, acceptor)
				end
			end
		}

		results.energy = energy
		return results
	end
end

