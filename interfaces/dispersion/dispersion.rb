################################################################################
#
# Interface Dispersion
#
# Author: Jan Rezac
# Date created: 2012-12-14
# License: Cuby4 license
# Description: Implementation of Jurecka's DFT-D dispersion
# Status: Works
#
################################################################################

#===============================================================================
# Implementation of Jurecka's DFT-D dispersion
#===============================================================================

require "classes/math/func"
require "classes/geometry/hybridization"
require "classes/calculation/default_method_settings"


class DispersionElement
	attr_accessor :c, :n, :r, :c_hyb, :r_hyb
end

module InterfaceDispersion
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient, :hessian]
	MODIFIER = true
	#=======================================================================

	#=======================================================================
	# Interface data - shared by all instances
	#=======================================================================
	@@data = nil

	#=======================================================================
	# Interface methods
	#=======================================================================

	def prepare_interface
		# Load data
		unless @@data
			# load the database
			Cuby::log.puts_debug "Dispersion: Loading element data"
			f = File.open(interface_dir + "/dispersion_data.yaml")
			@@data = YAML::load(f)["elements"]
			f.close

			# Load the parameters
			f = File.open(interface_dir + "/dispersion_parameters.yaml")
			@@parameters = YAML::load(f)
			f.close
		end

		# Determine parameters automatically from parent calculation setup
		default_settings = DefaultMethodSettings.new(@@parameters)
		method_defaults_loaded = default_settings.find_and_merge(@settings)

		# Check elements in system for having parameters
		@geometry.elements_in_system.each {|element|
			if @@data[element].c == 0.0
				Cuby::warning("No parameters found for element #{element.to_s} in system #{@name}. Dispersion correction will ignore these atoms.", :once)
			end
		}

		# Build atom r and c
		atomic_parameters

	end

	def calculate_interface
		# Load parameters from input
		parameters_from_input
		atomic_parameters_from_input

		# Calculation
		results = dispersion_energy_gradient(@what.include?(:gradient))

		# Hessian is calculated separately
		results.hessian = dispersion_hessian if @what.include?(:hessian)

		return results
	end

	#=======================================================================
	# Private methods: calculation
	#=======================================================================
	
	PM7_CUTOFF = 6.5

	def dispersion_energy_gradient(do_grad = false)
		# Initialize results
		results = Results.new

		# Damping function
		damp = Func::DampingFunction.new(1.0, @alpha, :df_0_to_1)

		# Development option - distance-dependent s6 for PM7
		if @settings[:dispersion_pm7_dd]
			if @settings.set?(:dispersion_pm7_r2)
				r1 = @settings[:dispersion_pm7_r1]
				r2 = @settings[:dispersion_pm7_r2]
			else
				r1 = @settings[:dispersion_pm7_r1] - @settings[:dispersion_pm7_w]
				r2 = @settings[:dispersion_pm7_r1] + @settings[:dispersion_pm7_w]
			end
			s6_long = @settings[:dispersion_pm7_s6d]
			dd_switch = Func::PolySwitching.new(r1,r2,:sw_1_to_0)
		end

		edisp = 0.0
		results.gradient = Gradient.zero(@geometry.size) if do_grad

		@geometry.each_index { |i|
			(i).times { |j|
				if ! (@geometry[i].properties[:ghost] || @geometry[j].properties[:ghost] || 
				      @geometry[i].properties[:dummy] || @geometry[j].properties[:dummy])

					if @settings[:dispersion_pm7_cutoff]
						next if @geometry[i].distance(@geometry[j]) >= PM7_CUTOFF
					end

					# calculate pairwise parameters
					r0,c6,params_ok = dispersion_get_pair_para(@geometry[i],@geometry[j], i, j)

					if params_ok
						# distance 
						rij = @geometry[i].distance(@geometry[j]) / 10.0	# in nanometers
						damp.r0 = @sr * r0
						eij = -c6/rij**6 * damp.calc(rij) / (1000.0*4.184*627.5095)

						if @settings[:dispersion_pm7_cutoff]
							eij *= 1.0 - Math::exp(-1.0*(rij*10 - PM7_CUTOFF)**2)
						end
						
						s6 = @s6
						if @settings[:dispersion_pm7_dd]
							s6 = s6_long + dd_switch.calc(rij*10.0) * (s6-s6_long)
						end

						edisp += eij * s6

						if do_grad
							# each coordinate: x,y,z
							3.times { |k|
								dcrd= (@geometry[i][k] - @geometry[j][k]) / 10.0 # dx in nanometers
								dgr = dcrd/rij * c6 * (-rij**-6 * damp.deriv(rij) + 6.0*rij**-7 * damp.calc(rij))
								dgr *= @s6
								# and apply to gradient
								results.gradient[i][k] += dgr / 1000.0 / 4.184 / NANOMETER2ANGSTROM 
								results.gradient[j][k] -= dgr / 1000.0 / 4.184 / NANOMETER2ANGSTROM
							}
						end
					end

				end
			}
		}

		results.energy = edisp * HARTREE2KCAL
		return results
	end

	#=======================================================================
	# Private methods: helpers
	#=======================================================================

	def parameters_from_input
		@alpha = @settings[:alpha]
		@sr = @settings[:sr]
		@s6 = @settings[:s6]
	end

	def atomic_parameters
		# Initialize parameter arrays
		@atom_r = []
		@atom_c = []

		# Evaluate hybridizations if needed
		if @settings[:dispersion_hyb]
			@hybridization = Hybridization.new(@geometry)
			@hybridization.evaluate_from_connectivity
		end

		@geometry.each_with_index{|atom, i|
			element = atom.element

			# Normal parameters
			c = @@data[element].c
			r = @@data[element].r

			if @settings[:dispersion_hyb]
				# Hybridization-dependent parameters
				if @@data[element].c_hyb
					c = @@data[element].c_hyb[@hybridization[i]]
				else
					Cuby::warning("Missing hybridization-dependent dispersion parameters for element #{element}", :once)
				end
				if @@data[element].r_hyb
					r = @@data[element].r_hyb[@hybridization[i]]
				else
					Cuby::warning("Missing hybridization-dependent dispersion parameters for element #{element}", :once)
				end
			end

			@atom_r[i] = r
			@atom_c[i] = c
		}
	end

	def atomic_parameters_from_input
		return unless @settings.set?(:dispersion_elements_r) ||
			      @settings.set?(:dispersion_elements_c)

		# Element-specific parameters from input
		user_r = @settings.elements_hash(:dispersion_elements_r)
		user_c = @settings.elements_hash(:dispersion_elements_c)

		@geometry.each_with_index{|atom, i|
			element = atom.element
			if @settings[:dispersion_hyb]
				# Hybridization-dependent parameters
				if r = user_r[element]
					@atom_r[i] = r[@hybridization[i]].to_f * 100.0 # Is in Angstrom
				end
				if c = user_c[element]
					@atom_c[i] = c[@hybridization[i]].to_f
				end
			else
				# Normal parameters, numbers expected
				if r = user_r[element]
					Cuby::error("Radius for element #{element} is not a number") unless r.kind_of?(Numeric)
					@atom_r[i] = r.to_f * 100.0 # Is in Angstrom
				end
				if c = user_c[element]
					Cuby::error("C6 coefficient for element #{element} is not a number") unless c.kind_of?(Numeric)
					@atom_c[i] = c.to_f
				end
			end
		}

	end

	def dispersion_get_pair_para(atom1,atom2, index1, index2)
		n1 = @@data[atom1.element].n
		n2 = @@data[atom2.element].n

		r1 = @atom_r[index1]
		r2 = @atom_r[index2]
		c1 = @atom_c[index1]
		c2 = @atom_c[index2]

		case @settings[:dispersion_mix] 
		when :grimme
			# Grimme: Arithmetic mean
			r0 = (r1 + r2) / 1000.0
			# C6
			c6 = 2.0 * c1 * c2 / (c1+c2)
		when :jurecka
			# cubic mean of Bondi vdW radii
			r0 = (r1**3 + r2**3) / (r1**2 + r2**2) / 1000.0 * 2.0	# in nanometers
			# Wu & Yang mixing of C6
			c6 = 2.0 * (c1**2 * c2**2 * n1 * n2)**(1.0/3.0)  /  ((c1 * n2**2)**(1.0/3.0) + (c2 * n1**2)**(1.0/3.0))
		when :jurecka_safe
			# cubic mean of Bondi vdW radii
			r0 = (r1**3 + r2**3) / (r1**2 + r2**2) / 1000.0 * 2.0	# in nanometers
			if  n1 == 0 || n2 == 0
				# Simple mixing when n's are missing
				c6 = 2.0 * c1 * c2 / (c1+c2)
			else
				# Wu & Yang mixing of C6
				c6 = 2.0 * (c1**2 * c2**2 * n1 * n2)**(1.0/3.0)  /  ((c1 * n2**2)**(1.0/3.0) + (c2 * n1**2)**(1.0/3.0))
			end
		else
			raise "Unknown mixing type \"#{@settings[:dispersion_mix].to_s}\""
		end
		c6 = 0.0 if c1 == 0 || c2 == 0
		ok = r1 > 0.0 && r2 > 0.0
		return [r0,c6,ok]
	end
end

__END__

	def dispersion_hessian
		hessian = Matrix.zero(@geometry.natom * 3)

		damp = Func::DampingFunction.new(1.0, @dispersion_para.alpha, :df_0_to_1)

		@geometry.each_index { |i|
			(i).times { |j|
			if ! (@geometry[i].properties[:ghost] || @geometry[j].properties[:ghost] || 
			      @geometry[i].properties[:dummy] || @geometry[j].properties[:dummy])

				# calculate pairwise parameters
				r0,c6,params_ok = dispersion_get_pair_para(@geometry[i],@geometry[j], i, j)

				if params_ok
					# distance 
					rij = @geometry[i].distance(@geometry[j]) / 10.0	# in nanometers

					damp.r0 = @dispersion_para.radii_scaling*r0
					deriv2 = c6*(-1.0/rij**6 * damp.deriv2(rij) +
						 	12.0/rij**7 * damp.deriv(rij) -
						 	42.0/rij**8 * damp.calc(rij))
					deriv1 = c6 * (-rij**-6 * damp.deriv(rij) + 6.0*rij**-7 * damp.calc(rij))
				
					3.times { |k|
					3.times { |l|
						a = (@geometry[i][k] - @geometry[j][k]) / 10.0 # dx in nanometers
						b = (@geometry[i][l] - @geometry[j][l]) / 10.0 # dx in nanometers

						fconst = deriv2 * a * b / rij**2 - deriv1 * a * b / rij**3
						fconst += deriv1 / rij if k == l


						fconst *= J2KCAL / NANOMETER2ANGSTROM / NANOMETER2ANGSTROM # convert from J/nm**2
						fconst *= @dispersion_para.global_scaling

						hessian[3*i+k, 3*j+l] -= fconst
						hessian[3*j+k, 3*i+l] -= fconst
						# fix i == j terms
						hessian[3*i+k, 3*i+l] += fconst 
						hessian[3*j+k, 3*j+l] += fconst 
					}
					}
				end
			end
			}
		}
		return hessian
	end
