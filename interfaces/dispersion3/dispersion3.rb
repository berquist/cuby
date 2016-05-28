################################################################################
#
# Interface Dispersion3
#
# Author: Jan Rezac
# Date created: 2012-12-14
# License: Cuby4 license
# Description: Implementation of Grimme's DFT-D3
# Status: Works
#
################################################################################

#===============================================================================
# Implementation of Grimme's DFT-D3
#===============================================================================

require "classes/geometry/hybridization"
require "classes/math/func"
require "classes/calculation/default_method_settings"
require "interfaces/dispersion3/classes/damping_t_t"

module InterfaceDispersion3
	#=======================================================================
	# Interface header
	#=======================================================================
	# Interface status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "OK"
	# Interface information
	INTERFACE = :calculation
	CAPABILITIES = [:energy, :gradient]
	IGNORE_AS_MODIFIER = [:atomic_charges]
	MODIFIER = true
	#=======================================================================

	#=======================================================================
	# C6 record class used in the YAML data file
	#=======================================================================
	class C6Record
		attr_accessor :c6
		attr_accessor :hyb_a
		attr_accessor :hyb_b
	end

	#=======================================================================
	# Interface data - shared by all intances
	#=======================================================================
	@@data = nil

	#=======================================================================
	# Interface methods
	#=======================================================================

	def prepare_interface
		# Load D3 data
		d3_data_load

		# Patch data
		d3_data_patch

		# Determine parameters automatically from parent calculation setup
		default_settings = DefaultMethodSettings.new(@@parameters)
		method_defaults_loaded = default_settings.find_and_merge(@settings)

		#: Check whether all parameters needed for a given calculation are set
		case @settings[:d3_damping]
		when :zero
			vars_to_check = [:d3_sr6, :d3_sr8, :d3_s6, :d3_s8, :d3_alpha6, :d3_alpha8]
		when :bj
			vars_to_check = [:d3_s6, :d3_s8, :d3_a1, :d3_a2]
		when :sapt
			vars_to_check = [:d3_s6, :d3_s8, :d3_a1, :d3_alpha6]
		when :tt
			vars_to_check = [:d3_s6, :d3_s8, :d3_a1, :d3_a2]
		when :none
			vars_to_check = [:d3_s6, :d3_s8]
		end
		vars_to_check << :d3_hh_para if @settings[:d3_hh_fix]
		missing = []
		vars_to_check.each{|varname|
			missing << varname.to_s unless @settings.set?(varname) || Settings.keyword(varname).when_not_present != :die
		}
		if missing.size > 0
			if method_defaults_loaded
				Cuby::error("Dispersion3 parameters missing although defaults were loaded for the present method:\n#{missing.join(', ')}")
			else
				Cuby::error("Dispersion3 parameters missing:\n#{missing.join(', ')}")
			end
		end
	end

	def calculate_interface
		# Calculate
		results = d3_dispersion_e_g(@what.include?(:gradient))
		return results
	end

	#=======================================================================
	# Private methods: Parameter setup
	#=======================================================================
	
	def d3_data_load
		# Load the data file if not available already
		unless @@data
			# Psych YAML parser is very slow and the data file is large
			# Therefore, we cache the serialized data separately
			if FileTest.exist?(Cuby::local_data_dir + "/dispersion3_data.dat")
				# Load the raw data
				Cuby::log.puts_debug "Dispersion3: cached data used"
				df = File.open(Cuby::local_data_dir + "/dispersion3_data.dat")
				@@data = Marshal.load(df)
				df.close
			else
				# load the database
				Cuby::log.puts_debug "Dispersion3: rebuilding data cache"
				f = File.open(interface_dir + "/dispersion3_data.yaml")
				@@data = YAML::load(f)
				f.close
				# dump the data to cache
				df = File.open(Cuby::local_data_dir + "/dispersion3_data.dat","w+")
				Marshal.dump(@@data, df)
				df.close
			end

			# Load the parameters
			f = File.open(interface_dir + "/dispersion3_parameters.yaml")
			@@parameters = YAML::load(f)
			f.close
		end
		
		# Make the data available locally
		@covalent_radii = @@data["covalent_radii"]
		@c6 = @@data["c6"]
		@r2r4 = @@data["r2r4"]
		@r0 = @@data["r0ab"]

		# Experimental: make all hybridization numbers integer
		#--------------------------------------------------
		# @c6.each_pair{|e1,h1|
		# 	h1.each_pair{|e2,arr|
		# 		arr.each_index{|i|
		# 			arr[i].hyb_a = arr[i].hyb_a.round.to_f
		# 			arr[i].hyb_b = arr[i].hyb_b.round.to_f
		# 		}
		# 	}
		# }
		#-------------------------------------------------- 
	end

	def d3_data_patch
		if @settings.set?(:d3_data_patch)
			f = File.open(@settings[:d3_data_patch])
			data_patch = YAML::load(f)
			f.close
			
			@c6_patch = data_patch["c6"]
		else
			@c6_patch = nil
		end
	end

	def d3_parameters_from_input
		#: Read parameters from input if present, override what has been
		#: set up previously

		@s6 = @settings[:d3_s6]
		@s8 = @settings[:d3_s8]

		case @settings[:d3_damping]
		when :zero
			@sr6 = @settings[:d3_sr6]
			@sr8 = @settings[:d3_sr8]
			@alpha = @settings[:d3_alpha6]
			if @settings.set?(:d3_alpha8)
				@alpha8 = @settings[:d3_alpha8]
			else
				@alpha8 = @alpha + 2.0
			end
		when :bj
			@a1 = @settings[:d3_a1]
			@a2 = @settings[:d3_a2]
		when :sapt
			@a1 = @settings[:d3_a1]
			@alpha = @settings[:d3_alpha6]
			if @settings.set?(:d3_alpha8)
				@alpha8 = @settings[:d3_alpha8]
			else
				@alpha8 = @alpha
			end
		when :tt
			@a1 = @settings[:d3_a1]
			@a2 = @settings[:d3_a2]
		end
		# H-H fix for SQM
		@hh_fix = @settings[:d3_hh_fix]
		if @hh_fix
			@hh_k =  @settings[:d3_hh_para][:k]
			@hh_e =  @settings[:d3_hh_para][:e]
			@hh_r =  @settings[:d3_hh_para][:r0]

			if @settings[:d3_hh_fix_version] == 2
				require "interfaces/dispersion3/classes/soft_core_repulsion"
				@hh_fix_2 = SoftCoreRepulsion.new(1.0,1.5, @settings[:d3_hh_para][:k], @settings[:d3_hh_para][:e], @settings[:d3_hh_para][:e2])
			end
		end
	end


	#=======================================================================
	# Private methods: Calculation
	#=======================================================================

	def d3_dispersion_e_g(grad = false)
	#: Dispersion energy calculation, complete -D3 form with continuous coordination
	
		# re-read parameters from input
		# (needed for parameterization where the input changes)
		d3_parameters_from_input 
		
		# Set up empty results
		results = Results.new
		results.gradient = Gradient.zero(@geometry.size) if grad

		# Hybridization
		hybridization = Hybridization.new(@geometry)
		case @settings[:d3_hybridization]
		when :fixed
			hybridization.evaluate_from_connectivity
		when :grimme
			if grad
				Cuby::recommendation("Gradient of the D3 correction with continuous hybridization not implemented.\nPlease set keyword 'd3_hybridization' to 'fixed'.")
			end
			# Using the radii from the data file
			hybridization.evaluate_grimme(@covalent_radii, false, grad)
			hybridization_grad = hybridization.gradient if grad
		when :grimme_mod
			if grad
				Cuby::recommendation("Gradient of the D3 correction with continuous hybridization not implemented.\nPlease set keyword 'd3_hybridization' to 'fixed'.")
			end
			# Using the radii from the data file
			hybridization.evaluate_grimme_mod(@covalent_radii, false, grad)
			hybridization_grad = hybridization.gradient if grad
		end

		# Cutoff
		if @settings[:d3_cutoff] != 0.0
			cutoff = true
			cutoff_r = @settings[:d3_cutoff]
			cutoff_roff = @settings[:d3_cutoff] + @settings[:d3_cutoff_buffer]
			cutoff_f = Func::PolySwitching.new(cutoff_r, cutoff_roff, :sw_1_to_0)
		else
			cutoff = false
		end

		# Store distances and c6 coefficients
		if @settings[:d3_3body]
			save_matrices = true
			matrix_r = Matrix.zero(@geometry.size)
			matrix_r0 = Matrix.zero(@geometry.size)
			matrix_c6 = Matrix.zero(@geometry.size)
		else
			save_matrices = false
		end

		# What to calculate?
		do_6 = @s6 != 0.0
		do_8 = @s8 != 0.0

		# What to print?
		print_pairwise = @settings[:d3_print].include?(:pairwise_energies)
		Cuby::log.puts "Pairwise D3 dispersion energies (kcal/mol):" if print_pairwise

		# Initialize sums
		edisp6 = 0.0
		edisp8 = 0.0
		erep = 0.0

		@geometry.each_index { |i|
			atom_a = @geometry[i]
			(i).times { |j|
				atom_b = @geometry[j]
				next if (atom_a.properties[:ghost] || atom_b.properties[:ghost] || 
					 atom_a.properties[:dummy] || atom_b.properties[:dummy])

				# Pairwise data
				r = atom_a.distance(atom_b)
				r0 = get_r0(atom_a.element, atom_b.element)
				if hybridization_grad
					# Gradients of hybridization number with respect to cartesian coordinates
					hg_i = hybridization_grad[i]
					hg_j = hybridization_grad[j]
					# Convert to derivatives with respect to r
					r_vec = (atom_a - atom_b)/r # Normal vector
					hg_i = hg_i.dot(r_vec)
					hg_j = -hg_j.dot(r_vec)
					c6, c6_deriv = interpolate_c6_with_deriv(atom_a.element, atom_b.element, hybridization[i], hybridization[j], hg_i, hg_j)
					c6_deriv = c6_deriv
				else
					c6 = interpolate_c6(atom_a.element, atom_b.element, hybridization[i], hybridization[j])
				end

				if save_matrices
					matrix_r[i,j] = r
					matrix_r0[i,j] = r0
					matrix_c6[i,j] = c6
				end

				# Default: zero
				e6 = 0.0
				e8 = 0.0
		
				# Cutoff		
				next if cutoff && r > cutoff_roff

				c8 = c8_from_c6(atom_a.element, atom_b.element, c6) 

				case @settings[:d3_damping]
				when :bj
					r0ab = (c8 / c6)**0.5
					fd = @a1*r0ab + @a2

					if do_6
						e6 = -@s6 * c6 / ((r * ANGSTROM2BOHR)**6 + fd**6) * HARTREE2KCAL
						if grad
							deriv_6 = 6.0 * c6 * r**5 * @s6 * ANGSTROM2BOHR**6 * HARTREE2KCAL / (r**6 * ANGSTROM2BOHR**6 + fd**6)**2
							deriv = deriv_6
						end
					end

					if do_8
						e8 = -@s8 * c8 / ((r * ANGSTROM2BOHR)**8 + fd**8) * HARTREE2KCAL
						if grad
							deriv_8 = 8.0 * c8 * r**7 * @s8 * ANGSTROM2BOHR**8 * HARTREE2KCAL / (r**8 * ANGSTROM2BOHR**8 + fd**8)**2
							deriv += deriv_8
						end
					end
				when :zero
					if do_6
						damping6 = damping_function(r, r0, @sr6, @alpha)
						e6 = -@s6 * c6 / (r * ANGSTROM2BOHR)**6 *  HARTREE2KCAL
						if grad
							deriv_damp6 = damping_function_deriv(r, r0, @sr6, @alpha)

							if hybridization_grad
								e6_deriv = -@s6 * HARTREE2KCAL / ANGSTROM2BOHR**6 *
									(c6 * -6.0/(r)**7 + c6_deriv * 1.0/r**6)
							else
								e6_deriv = -@s6 * c6 * HARTREE2KCAL * -6.0/(r)**7 / ANGSTROM2BOHR**6
							end
							deriv6 = e6 * deriv_damp6 + e6_deriv * damping6
							deriv = deriv6

						end
						e6 *= damping6
					end

					if do_8
						e8 = -@s8 * c8 / (r * ANGSTROM2BOHR)**8 *  HARTREE2KCAL
						damping8 = damping_function(r, r0, @sr8, @alpha8)

						if grad
							deriv_damp8 = damping_function_deriv(r, r0, @sr8, @alpha8)
							e8_deriv = -@s8 * c8 * HARTREE2KCAL * -8.0/(r)**9 / ANGSTROM2BOHR**8
							deriv8 = e8 * deriv_damp8 + e8_deriv * damping8
							deriv += deriv8
						end
						e8 *= damping8
					end

				when :sapt
					damping = Math::erf(@a1 * r / r0)**@alpha
					damping8 = Math::erf(@a1 * r / r0)**@alpha8
					# **3, a1 = 1.19
					# **4, a1 = 1.3, RMSE CT dataset 0.800
					# **5, a1 = 1.37, RMSE 0.679
					# **6, a1 = 1.44, RMSE 0.669

					e6 = -@s6 * c6 / (r * ANGSTROM2BOHR)**6 * HARTREE2KCAL
					e8 = -@s8 * c8 / (r * ANGSTROM2BOHR)**8 * HARTREE2KCAL

					if grad
						Cuby::error "Gradient not implemented for SAPT damping"
					end

					e6 *= damping
					e8 *= damping8
				when :tt
					#r0ab = (c8 / c6)**0.5
					r0ab = r0

					e6 = -@s6 * c6 / (r * ANGSTROM2BOHR)**6 *  HARTREE2KCAL
					e8 = -@s8 * c8 / (r * ANGSTROM2BOHR)**8 *  HARTREE2KCAL

					if grad
						Cuby::error "Gradient not implemented for Tang and Toennies damping"
					end

					e6 *= DampingTT.f6(r,r0ab * @a1 + @a2)
					e8 *= DampingTT.f8(r,r0ab * @a1 + @a2)
				when :none
					e6 = -@s6 * c6 / (r * ANGSTROM2BOHR)**6 *  HARTREE2KCAL if do_6
					e8 = -@s8 * c8 / (r * ANGSTROM2BOHR)**8 *  HARTREE2KCAL if do_8

					if grad
						e6_deriv = -@s6 * c6 * HARTREE2KCAL * -6.0/(r)**7 / ANGSTROM2BOHR**6
						deriv = e6_deriv

						e8_deriv = -@s8 * c8 * HARTREE2KCAL * -8.0/(r)**9 / ANGSTROM2BOHR**8
						deriv += e8_deriv
					end
				end

				# Hydrogen - hydrogen fix for SQM
				if @hh_fix && atom_a.element == :H && atom_b.element == :H
					if @hh_fix_2
						erep += @hh_fix_2.calculate(r)
						deriv += @hh_fix_2.deriv(r) if grad
					else
						# energy
						erep += @hh_k * (1.0 - 1.0/(1.0 + Math::exp(-@hh_e*(r/@hh_r - 1.0))))
						# gradient
						if grad
							erep_d = (1.0 / (1.0 + Math::exp(-@hh_e*(r/@hh_r-1.0)))**2 * @hh_e/@hh_r * Math::exp(-@hh_e*(r/@hh_r-1.0))) * -@hh_k
							deriv += erep_d
						end
					end
				end

				# apply cutoff
				if cutoff &&  r > cutoff_r
					cut = cutoff_f.calc(r)
					if grad
						cut_d = cutoff_f.deriv(r) 
						deriv = deriv * cut + cut_d * (e6+e8)
					end
					e6 *= cut
					e8 *= cut
				end

				if print_pairwise
					Cuby::log.puts(sprintf("%8d%8d%12.6f",i+1,j+1,e6+e8))
				end

				edisp6 += e6
				edisp8 += e8

				if grad
					# Project gradient size to coordinate
					dcrd= @geometry[i] - @geometry[j]
					dgr = deriv/r * dcrd

					# and apply to gradient
					results.gradient[i] += dgr
					results.gradient[j] -= dgr
				end
			}
		}

		results.energy_components[:dispersion_c6] = edisp6
		results.energy_components[:dispersion_c8] = edisp8
		results.energy_components[:hh_repulsion] = erep if @hh_fix

		edisp =  edisp6 + edisp8 + erep

		# 3-body term, calculated separately
		if @settings[:d3_3body]
			if grad
				Cuby::error("Gradient of the D3 correction with 3-body term not implemented.")
			end
			if cutoff
				Cuby::error("Cutoff in D3 correction with 3-body term not implemented.")
			end
			edisp3b = d3_3body_e(hybridization, matrix_r, matrix_r0, matrix_c6) * @settings[:d3_3body_scaling]
			results.energy_components[:dispersion_3body] = edisp3b
			edisp += edisp3b
		end

		results.energy = edisp * @settings[:d3_scaling]
		return results
	end

	def d3_3body_e(hybridization, matrix_r, matrix_r0, matrix_c6)
		# Settings
		@sr3b = @settings[:d3_3body_sr]

		edisp3b = 0.0
		@geometry.each_index { |i| 
			next if (@geometry[i].properties[:ghost] || @geometry[i].properties[:dummy])
			atom_a = @geometry[i]
			(i).times { |j| 
				next if (@geometry[j].properties[:ghost] || @geometry[j].properties[:dummy])
				atom_b = @geometry[j]
				(j).times { |k|
					next if (@geometry[k].properties[:ghost] || @geometry[k].properties[:dummy])
					atom_c = @geometry[k]

					# Distances
					rab = matrix_r[i,j]
					rbc = matrix_r[j,k]
					rca = matrix_r[i,k]

					# Angles
					aa = Coordinate.angle(atom_c, atom_a, atom_b)
					ab = Coordinate.angle(atom_a, atom_b, atom_c)
					ac = Coordinate.angle(atom_b, atom_c, atom_a)

					# C9 coef
					c9 = -1.0 * (matrix_c6[i,j] * matrix_c6[j,k] * matrix_c6[i,k])**0.5

					# r0
					r0ab = matrix_r0[i,j]
					r0bc = matrix_r0[j,k]
					r0ca = matrix_r0[i,k]

					# Damping function
					r_avg = (rab * rbc * rca)**(1.0/3.0)
					r0_avg = (r0ab * r0bc * r0ca)**(1.0/3.0)
				
					if @sr3b == 0
						# no damping
						fd = 1.0
					else	
						fd = damping_function(r_avg, r0_avg, @sr3b, 16.0)
					end

					if @settings[:development][:damping_3b_tt]
						a = @settings[:x0]
						b = @settings[:x1]
						fd =  DampingTT.f6(rab,r0ab * a + b) *
						      DampingTT.f6(rbc,r0bc * a + b) *
						      DampingTT.f6(rca,r0ca * a + b)
					elsif @settings[:development][:damping_3b_erf]
						a = @settings[:x0]
						alpha = @settings[:x1]
						fd = Math::erf(a * rab / r0ab)**alpha *
						     Math::erf(a * rbc / r0bc)**alpha *
						     Math::erf(a * rca / r0ca)**alpha
					end

					# Dispersion energy
					edisp3b -= fd * c9 * (3.0 * Math::cos(aa) * Math::cos(ab) * Math::cos(ac) + 1.0) / (rab*rbc*rca)**3 / (ANGSTROM2BOHR)**9 * HARTREE2KCAL
				}
			}
		}
		return edisp3b
	end

	#=======================================================================
	# Private methods: misc
	#=======================================================================
	
	def interpolate_c6(ele_a, ele_b, cn_a, cn_b)
		z = 0.0
		w = 0.0

		if @c6_patch && @c6_patch[ele_a] && @c6_patch[ele_a][ele_b]
			Cuby::log.puts_debug("c6 for #{ele_a} - #{ele_b} taken from patch file")
			@c6_patch[ele_a][ele_b].each{|c6rec|
				l = Math::exp(-1.0 * 4.0 * ((cn_a - c6rec.hyb_a)**2 + (cn_b - c6rec.hyb_b)**2))
				f = 1.0
				z += c6rec.c6 * f * l
				w += l
			}
		elsif @c6_patch && @c6_patch[ele_b] && @c6_patch[ele_b][ele_a]
			return interpolate_c6(ele_b, ele_a, cn_b, cn_a)
		elsif @c6[ele_a][ele_b]
			@c6[ele_a][ele_b].each{|c6rec|
				l = Math::exp(-4.0 * ((cn_a - c6rec.hyb_a)**2 + (cn_b - c6rec.hyb_b)**2))
				z += c6rec.c6 * l
				w += l
				unless @settings[:d3_legacy_version]
					# Same element - symmetric points are missing in the table, must be constructed here
					if ele_a == ele_b && c6rec.hyb_a != c6rec.hyb_b
						l = Math::exp(-4.0 * ((cn_a - c6rec.hyb_b)**2 + (cn_b - c6rec.hyb_a)**2))
						z += c6rec.c6 * l
						w += l
					end
				end
			}
		elsif @c6[ele_b][ele_a]
			return interpolate_c6(ele_b, ele_a, cn_b, cn_a)
		else
			raise "No C6 parameters found for pair #{ele_a} ... #{ele_b}"
		end

		return z/w
	end

	def interpolate_c6_with_deriv(ele_a, ele_b, cn_a, cn_b, dcn_a, dcn_b)
		z = 0.0
		w = 0.0
		dz = 0.0
		dw = 0.0

		if @c6[ele_a][ele_b]
			@c6[ele_a][ele_b].each{|c6rec|
				l = Math::exp(-1.0 * 4.0 * ((cn_a - c6rec.hyb_a)**2 + (cn_b - c6rec.hyb_b)**2))
				dl = -4.0 * l * (2.0*(cn_b - c6rec.hyb_b)*dcn_b + 2.0 *(cn_a - c6rec.hyb_a)*dcn_a)
				z += c6rec.c6 * l
				w += l
				dz += c6rec.c6 * dl
				dw += dl
				unless @settings[:d3_legacy_version]
					# Same element - symmetric points are missing in the table, must be constructed here
					if ele_a == ele_b && c6rec.hyb_a != c6rec.hyb_b
						l = Math::exp(-1.0 * 4.0 * ((cn_a - c6rec.hyb_b)**2 + (cn_b - c6rec.hyb_a)**2))
						dl = -4.0 * l * (2.0*(cn_b - c6rec.hyb_a)*dcn_b + 2.0 *(cn_a - c6rec.hyb_b)*dcn_a)
						z += c6rec.c6 * l
						w += l
						dz += c6rec.c6 * dl
						dw += dl
					end
				end
			}
		elsif @c6[ele_b][ele_a]
			c6, deriv = interpolate_c6_with_deriv(ele_b, ele_a, cn_b, cn_a, dcn_b, dcn_a)
			return [c6, deriv]
		else
			raise "No C6 parameters found for pair #{ele_a} ... #{ele_b}"
		end

		derivative = (dz*w - z*dw)/w**2

		return [z/w, derivative]
	end

	def c8_from_c6(ele_a, ele_b, c6)
		s42 = 0.5 # Constant

		q_a = s42 * PeriodicTable.proton_number(ele_a).to_f**0.5 * @r2r4[ele_a]
		q_b = s42 * PeriodicTable.proton_number(ele_b).to_f**0.5 * @r2r4[ele_b]

		c8 = 3.0 * c6 * (q_a * q_b)**0.5
		return c8
	end

	def damping_function(r, r0, sr, alpha)
		return 1.0 / (1.0 + 6.0 * (r/(sr*r0))**(-1.0 * alpha))
	end

	def damping_function_deriv(r, r0, sr, alpha)
		z = (r/(sr*r0))**alpha
		return alpha * z / (r *(z+6.0)) - alpha * z**2/ (r*(z+6.0)**2)
	end

	def get_r0(ele_a, ele_b)
		if @r0[ele_a][ele_b]
			return @r0[ele_a][ele_b]
		else
			return @r0[ele_b][ele_a]
		end
	end
end
