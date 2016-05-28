################################################################################
#
# Module HBonds4
#
# Author: J. Rezac
# Date created: 2010-12-14 
# License: Cuby license
# Description: Fourth generation H-bonding correction
# Extends: Calculation
# Status: Not documented, in progress
#
################################################################################

require "classes/math/polynomial_from_conditions"
require "classes/math/func"
require "classes/geometry/continuous_valence"
require "classes/calculation/default_method_settings"

module InterfaceHBonds4
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "All cuby3 functionality implemented"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient]
	MODIFIER = true
	#=======================================================================

	#=======================================================================
	# Interface methods
	#=======================================================================

	def prepare_interface
		h_bonds_load_parameters
	end

	def calculate_interface
		# Calculate
		results = h_bonds_energy_grad(@what.include?(:gradient))
		return results
	end

	#=======================================================================
	# Constants
	#=======================================================================
	# Possible donors and acceptors (Cuby selection)
	HB_DONOR_ACCEPTOR = "@N,O"
	# Distances determining the potential curve
	HB_R_0 = 1.5
	HB_R_MIN = 3.0
	HB_R_CUTOFF = 5.5
	# Angular cutoff
	HB_A_CUTOFF = Math::PI / 2

	@@hb_data = nil

	#=======================================================================
	# Private methods
	#=======================================================================
	
	def h_bonds_load_parameters
		# Build list of parameters
		# First index is donor, second acceptor
		@parameters = {
			:O => {},
			:N => {},
		}

		# Load default parameters from yaml file
		unless @settings.set?(:h_bonds4_parameters)
			unless @@hb_data
				@@hb_data = YAML::load(File.open(interface_dir + "/parameters.yaml"))
			end

			# Method lookup
			default_settings = DefaultMethodSettings.new(@@hb_data)
			unless default_settings.find_and_merge(@settings)
				Cuby::error("Default parameters for H-bonds4 not found for this method")
			end
		end

		# Load parameters from settings
		@parameters[:O][:O] = @settings[:h_bonds4_parameters]["oh_o"]
		@parameters[:O][:N] = @settings[:h_bonds4_parameters]["oh_n"]
		@parameters[:N][:O] = @settings[:h_bonds4_parameters]["nh_o"]
		@parameters[:N][:N] = @settings[:h_bonds4_parameters]["nh_n"]

		# Generate the polynomial for the radial term from a set of conditions
		@hb_func = PolynomialFromConditions.build([
			PolynomialFromConditions::Condition.new(:value, HB_R_0, 0.0),
			PolynomialFromConditions::Condition.new(:value, HB_R_MIN, -1.0),
			PolynomialFromConditions::Condition.new(:value, HB_R_CUTOFF, 0.0),
			PolynomialFromConditions::Condition.new(:first_deriv, HB_R_0),
			PolynomialFromConditions::Condition.new(:first_deriv, HB_R_MIN),
			PolynomialFromConditions::Condition.new(:first_deriv, HB_R_CUTOFF),
			PolynomialFromConditions::Condition.new(:second_deriv, HB_R_CUTOFF),
			PolynomialFromConditions::Condition.new(:third_deriv, HB_R_CUTOFF)
		])
		# Save the derivative of the radial term
		@hb_func_deriv = @hb_func.derivative

		# Make sure the polynomial has only one extreme (minimum) in the range used - Development only
		#--------------------------------------------------
		# extremes = []
		# @hb_func_deriv.roots.each{|u|
		# 	if u.imag.abs < 1.0e-6 # Allow for some numerical error
		# 		if u.real > HB_R_0 + 1.0e-4 && u.real < HB_R_CUTOFF - 1.0e-4
		# 			extremes << u
		# 		end
		# 	end
		# }
		# if extremes.count > 1
		# 	message = "H-bond correction potential has more than one extreme:\n#{extremes.map{|x| x.to_s}.join("\n")}"
		# 	Cuby::error(message)
		# end
		#-------------------------------------------------- 

		# Prepare list of donors and acceptors
		@d_a_list = @geometry.atomlist_from_selection(HB_DONOR_ACCEPTOR)
		@h_list = @geometry.atomlist_from_selection("@H")
	end

	def h_bonds_energy_grad(grad = false)
		# Initialize energy and gradient wit zeroes
		results = Results.new
		hb_energy = 0.0
		results.gradient = Gradient.zero(@geometry.size) if grad

		# Evaluate valence for all atoms
		continuous_valence = ContinuousValence.new(@geometry, grad, true)

		# Set up angular function
		@f_ang = Func::PolySwitching.new(0, HB_A_CUTOFF)
		
		# Iterate over all possible donor/acceptor pairs
		@d_a_list.each_index{|ii|
			ii.times{|jj|
				atom_a = @geometry[@d_a_list[ii]]
				atom_b = @geometry[@d_a_list[jj]]

				# Calculate distance, skip cycle if distance outside interval where radial term is nonzero
				rab = atom_a.distance(atom_b)
				next if rab > HB_R_CUTOFF
				next if rab < 1.5

				# Iterate over all hydrogens
				@h_list.each{|hi|
					atom_h = @geometry[hi]
					rah = atom_h.distance(atom_a)
					rbh = atom_h.distance(atom_b)

					# Pythagorean theorem used to pick hydrogens with angle XHY > 90 deg
					# This also filter distant hydrogens
					next if rah**2 + rbh**2 >= rab**2

					# Donor and acceptor definition
					# Donor is the atom closer to the hydrogen
					h = atom_h
					if atom_a.distance(h) < atom_b.distance(h)
						index_d = @d_a_list[ii]
						index_a = @d_a_list[jj]
						donor = atom_a
						acceptor = atom_b
					else
						index_d = @d_a_list[jj]
						index_a = @d_a_list[ii]
						donor = atom_b
						acceptor = atom_a
					end

					# Hack for skipping oxygen acceptors in O=S groups
					if acceptor.element == :O && @settings[:h_bonds4_skip_acceptor].include?("OS")
						sulphur = false
						continuous_valence.atoms[index_a].each{|i|
							atom = @geometry.at(i)
							if atom.element == :S && atom.distance(acceptor) < 1.1 * (PeriodicTable.covalent_radius(:S)+PeriodicTable.covalent_radius(:O))
								sulphur = true
								break
							end
						}
						next if sulphur
					end

					# Bond switching function
					rdh = donor.distance(h)
					rah = acceptor.distance(h)

					if rdh > 1.15
						# Hydrogen is leaving - correction is scaled by the bond switching function
						rdhs = rdh - 1.15
						ravgs = 0.5*rdh + 0.5*rah - 1.15
						x = rdhs/ravgs

						#--------------------------------------------------
						# bsfunc = Func::PolySwitching.new(0.0, 1.0, :sw_1_to_0)
						# bond_switch = bsfunc.calc(x)
						# bond_switch_d = bsfunc.deriv(x)
						#-------------------------------------------------- 
						# Rewritten explicitly as:
						bond_switch = 1.0-(-20.0*x**7 + 70.0*x**6 -84.0*x**5 + 35.0*x**4)
						bond_switch_d = -(-140.0*x**6 + 420.0*x**5 - 420.0*x**4 + 140.0*x**3)

						if grad
							# Bond switching derivatives
							dhd = bond_switch_d / ravgs
							dhda = 0.5 * bond_switch_d * -x / ravgs

							dbs_d = (donor-atom_h)/rdh * dhd + (donor-atom_h)/rdh * dhda
							dbs_a = (acceptor-atom_h)/rah * dhda
							dbs_h = -dbs_d + -dbs_a
						end
					else
						# Hydrogen is within covent distance from donor, no scaling
						bond_switch = 1.0
						if grad
							dbs_h = dbs_d = dbs_a = Coordinate.new
						end
					end

					# Radial function value
					poly = @hb_func.calc(rab)
					# ... and derivatives
					if grad
						dpoly = @hb_func_deriv.calc(rab)
						drad_d = (donor - acceptor)/rab * dpoly
						drad_a = -drad_d
					end

					# Angular
					angle = Math::PI - Coordinate.angle(donor, h, acceptor) # Angle is defined so that it is 0 in linear XHY arrangement
					next if angle >= HB_A_CUTOFF
					angular = 1.0 - @f_ang.calc(angle)**2
					# ... and derivatives
					if grad
						if angle == 0.0
							dang_d = dang_a = dang_h = Coordinate.new
						else
							dangular = -@f_ang.deriv(angle) * 2.0 * @f_ang.calc(angle)
							u = donor - h
							v = acceptor - h
							dot = u.dot(v)
							dang_d = -dangular / (1.0 - dot**2 / rdh**2 / rah**2)**0.5 * -(v/rdh/rah - u*dot/rdh**3/rah)
							dang_a = -dangular / (1.0 - dot**2 / rdh**2 / rah**2)**0.5 * -(u/rdh/rah - v*dot/rdh/rah**3)
							dang_h = -dang_d - dang_a
						end
					end
					
					# Lookup parameters
					c = @parameters[donor.element][acceptor.element]


					# Multiplying factor and gradients for scaling specific atom types
					scaling_mul = 1.0
					scaling_grads = {} # Hash of hashes, first index is central atom, second index is each atom paired with the central one
					scaling_factors = {} ## first index is central atom, second is the associated scaling factor

					# Individual cases of scaling
					
					if @settings[:h_bonds4_scale_charged]
						# COO- acceptor
						if acceptor.element == :O
					 		mult_slope =  @settings[:h_bonds4_parameters]["multiplier_coo"] - 1.0

							o1 = acceptor
							cc = nil
							o2 = nil

							o1_i = index_a
							cc_i = nil
							o2_i = nil

							# Look for c and o2
							cdist = 9.9e9
							continuous_valence.atoms[o1_i].each{|i|
								atom = @geometry[i]
								if atom.element == :C
									dist = o1.distance(atom)
									if dist < cdist
										cc = atom
										cc_i = i
										cdist = dist
									end
								end
							}
							if cc # We have C atom close
								odist = 9.9e9
								continuous_valence.atoms[cc_i].each{|i|
									atom = @geometry[i]
									if atom.element == :O && atom != o1
										dist = cc.distance(atom)
										if dist < odist
											o2 = atom
											o2_i = i
											odist = dist
										end
									end

								}
							end
							if o2 # The three atoms in COO
								# o1
								v = continuous_valence.valence[o1_i]
								v = (1.0 - (1.0 - v).abs)
								v = 0.0 if v < 0.0
								vo1 = v

								# o2
								v = continuous_valence.valence[o2_i]
								v = (1.0 - (1.0 - v).abs)
								v = 0.0 if v < 0.0
								vo2 = v

								# cc
								v = continuous_valence.valence[cc_i]
								v = (1.0 - (3.0 - v).abs)
								v = 0.0 if v < 0.0
								vcc = v

								scaling_mul *= 1.0 + mult_slope * vo1 * vo2 * vcc

								if grad
									# O1
									if vo1 != 0.0
										grad_hash_o1 = {}
										continuous_valence.gradient[o1_i].each_pair{|i,g|
											atom = @geometry[i]
											x = o1.distance(atom)
											g *= -1.0 if continuous_valence.valence[o1_i] > 1.0
											grad_hash_o1[i] = (o1 - atom) * -g/x * mult_slope * vo2 * vcc
										}
										scaling_grads[o1_i] = grad_hash_o1
										scaling_factors[o1_i] = 1.0 + mult_slope * vo2 * vcc * vo1
									end
									
									# O2
									if vo2 != 0.0
										grad_hash_o2 = {}
										continuous_valence.gradient[o2_i].each_pair{|i,g|
											atom = @geometry[i]
											x = o2.distance(atom) 
											g *= -1.0 if continuous_valence.valence[o2_i] > 1.0
											grad_hash_o2[i] = (o2 - atom) * -g/x * mult_slope * vo1 * vcc
										}
										scaling_grads[o2_i] = grad_hash_o2
										scaling_factors[o2_i] = 1.0 + mult_slope * vo1 * vcc * vo2
									end

									# CC
									if vcc != 0.0
										grad_hash_cc = {}
										continuous_valence.gradient[cc_i].each_pair{|i,g|
											atom = @geometry[i]
											x = cc.distance(atom)
											g *= -1.0 if continuous_valence.valence[cc_i] > 3.0
											grad_hash_cc[i] = (cc - atom) * -g/x * mult_slope * vo1 * vo2
										}
										scaling_grads[cc_i] = grad_hash_cc
										scaling_factors[cc_i] = 1.0 + mult_slope * vo1 * vo2 * vcc
									end
								end

							end
							
						end
					end

					if @settings[:h_bonds4_scale_charged]
						# NHR3+ donor
						# Tested
						if donor.element == :N
							mult_slope = @settings[:h_bonds4_parameters]["multiplier_nh4"] - 1.0
							hyb =  continuous_valence.valence[index_d]

							if hyb >= 3.0
								hyb -= 3.0
							else
								hyb = 0.0
							end

							if (grad && hyb != 0)
								grad_hash = {}
								continuous_valence.gradient[index_d].each_pair{|i,g|
									atom = @geometry[i]
									x = donor.distance(atom)
									grad_hash[i] = (donor - atom) * -g/x * mult_slope
								}
								scaling_grads[index_d] = grad_hash
								scaling_factors[index_d] = 1.0 + mult_slope * hyb
							end

							scaling_mul *= 1.0 + mult_slope * hyb
						end
					end

					# Water donor / O acceptor
					# Tested on proton transfer in water dimer
					if donor.element == :O && acceptor.element == :O
						# Count hydrogens and other atoms in vicinity
						hydrogens = 0.0
						others = 0.0
						continuous_valence.atoms[index_d].each{|i|
							atom = @geometry[i]
							if atom.element == :H
								hydrogens += continuous_valence.contribution(donor, atom)
							else
								others += continuous_valence.contribution(donor, atom)
							end
						}

						# Water-like atom
						if hydrogens > 1.0
							mult_slope = @settings[:h_bonds4_parameters]["multiplier_wh_o"] - 1.0
							v = hydrogens

							# Increase when number of hydrogens gets to 2
							if v > 1.0 && v <= 2.0
								v = v - 1.0
								sign = 1.0
							elsif v > 2.0 && v < 3.0
								v = 3.0 - v
								sign = -1.0
							else
								v = 0.0
							end

							# Reduce if other atoms are present
							v2 = 1.0 - others
							v2 = 0.0 if v2 < 0.0

							scaling_mul *= 1.0 + mult_slope * v * v2

							if grad && v != 0 && v2 != 0
								grad_hash = {}
								continuous_valence.atoms[index_d].each{|i|
									atom = @geometry[i]
									x = donor.distance(atom)
									if atom.element == :H
										g = continuous_valence.contribution_d(donor, atom)
										grad_hash[i] = (donor - atom) * -g/x * mult_slope * sign
									else
										g = continuous_valence.contribution_d(donor, atom)
										grad_hash[i] = (donor - atom) * g/x * mult_slope
									end
								}
								scaling_grads[index_d] = grad_hash
								scaling_factors[index_d] = 1.0 + mult_slope * v * v2
							end
						end
					end

					# Additional scaling
					scaling2 = 1.0
					if @settings.set?(:h_bonds4_extra_scaling)
						@settings[:h_bonds4_extra_scaling].each_pair{|selection, scaling_f|
							list = @geometry.atomlist_from_selection(selection)
							if list.include?(index_a)
								scaling2 *= scaling_f
							end
							if list.include?(index_d)
								scaling2 *= scaling_f
							end
						}
					end
					scaling_mul *= scaling2


					# Calculate energy
					e = c * poly * angular * bond_switch * scaling_mul
					# Add it to the total correction energy
					hb_energy += e

					# Derivatives
					if grad
						grad_d = dbs_d * poly * angular * c * scaling_mul
						grad_h = dbs_h * poly * angular * c * scaling_mul
						grad_a = dbs_a * poly * angular * c * scaling_mul

						grad_d += drad_d * angular * bond_switch * c * scaling_mul
						grad_a += drad_a * angular * bond_switch * c * scaling_mul

						grad_d += dang_d * poly * bond_switch * c * scaling_mul
						grad_h += dang_h * poly * bond_switch * c * scaling_mul
						grad_a += dang_a * poly * bond_switch * c * scaling_mul

						# Apply gradients from atomtype scaling
						scaling_grads.each_pair{|center,list|
								list.each_pair{|i,g|
									gc = g * (poly * angular * bond_switch * c)
									gc *= scaling_mul / scaling_factors[center] # multiply gradient with all but its own scaling factors
									results.gradient[i] += gc
									results.gradient[center] -= gc
								}
						}

						# Add this H-bond to the total gradient
						results.gradient[index_d] += grad_d * scaling2
						results.gradient[hi] += grad_h * scaling2
						results.gradient[index_a] += grad_a * scaling2
					end

					# Printing
					if false
						puts "Atoms: #{index_d + 1} #{hi + 1} #{index_a + 1}"
						puts "radial: #{poly}"
						puts "angular: #{angular}"
						puts "bond_switch: #{bond_switch}"
						puts "scaling: #{scaling_mul}"
						puts "energy: #{e}"
						puts "--------"
					end
				}

			}
		}

		if @settings[:h_bonds4_pt_corr] != 0.0
			f = @settings[:h_bonds4_pt_corr]
			# Iterate over all hydrogens
			@h_list.each{|hi|
				atom_h = @geometry[hi]
				cv = continuous_valence.valence[hi]
				cv = 1.0 if cv >= 1.0
				hb_energy += (1.0 - cv) * f

				# gradient
				if cv < 1.0 && grad
					continuous_valence.atoms[hi].each{|i|
						atom = @geometry[i]
						x = atom_h.distance(atom)
						g = continuous_valence.contribution_d(atom_h, atom) * f
						grad = (atom_h - atom) * g/x
						results.gradient[i] += grad
						results.gradient[hi] -= grad
					}
				end
			}
		end

		results.energy = hb_energy
		return results
	end
end

