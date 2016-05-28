require "classes/objects/vibrations.rb"
require "classes/objects/thermodynamics.rb"
require "classes/misc/array_of_arrays.rb"

module ProtocolFrequencies
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :elementary_calculation
	#=======================================================================
	
	attr_reader :vibrations
	attr_reader :thermodynamics

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Harmonic vibrations calculation')
	end

	def prepare_protocol
		# Read geometry
		load_geometry

		if @settings[:hessian_read] == ""
			# Create calculation
			@calculation = Calculation.new(@name, @settings, @geometry)
			# Prepare calculation
			@calculation.prepare([:energy, :gradient, :hessian])
		end
	end

	def run_protocol(queue)
		# Submit calculation
		if @settings[:hessian_read] == ""
			@calculation.send_to_queue(queue)
		end
	end

	def print_protocol
		if @settings[:hessian_read] == ""
			# Print results - energy
			results.print_energy
			if @settings[:print].include?(:energy_decomposition)
				Cuby::log.puts_v_only(:normal,"")
				results.print_energy_decomposition 
			end
			# Print results - gradient RMS
			Cuby::log.puts_v_only(:normal,"")
			results.gradient.project_out_trans_rot!(@geometry)
			Cuby::log.puts_v(:brief, "RMS grad: #{'%.6f' % results.gradient.rms} kcal/mol/A")
			# Save the hessian to a YAML file
			unless @settings[:hessian_write] == ""
				File.open(@settings[:hessian_write],"w+"){|f|
					f.puts results.hessian.to_yaml
				}
			end
		else
			# Method 'results' must be called before @vibrations and @thermodynamics are used
			results
		end


		# Print vibrational frequencies and related data
		Cuby::log.puts_v_only(:normal,"")
		@vibrations.print(@settings)

		# Print thermodynamic variables (in ideal gas/rigid rotor/harmonic oscillator approximation)
		Cuby::log.puts("")
		@thermodynamics.print(0.0) # Set electronic enerrgy to 0

		# Write molden input
		if @settings[:freq_molden_file] != ""
			Cuby::log.puts_debug("Writing vibrational data for molden to file '#{@settings[:freq_molden_file]}'")
			@vibrations.write_molden_input(@geometry, @settings, @settings[:freq_molden_file])
		end

		if @settings[:freq_print].include?(:transformation)
			Cuby::log.puts("")
			Cuby::log.puts("Transformation matrix - dimensionless q -> cartesian diplacement in A")
			Cuby::log.puts(transformation_normal_to_cart(@vibrations).to_s("%20.15e", " ", "", ""))
		end

		if @settings[:freq_mode_scan]
			if @settings[:freq_mode_scan_multidimensional].size == 0
				one_dimensional_scans(@vibrations)
			else
				multi_dimensional_scan(@vibrations, @settings[:freq_mode_scan_multidimensional])
			end
		end

	end

	def cleanup_protocol
		if @settings[:hessian_read] == ""
			@calculation.cleanup
		end
	end

	def results
		# Results are collected only first time they are needed
		unless @results
			# Get results
			if @settings[:hessian_read] == ""
				@results = @calculation.results
			else
				# When Hessian is read, set up empty results
				@results = Results.new
			end

			if @settings[:hessian_read] == ""
				# When the hessian is from a calculation, look also for dipole derivatives
				if results.multipoles && results.multipoles[:dipole] && results.multipoles[:dipole].derivative
					dipole_derivatives = results.multipoles[:dipole].derivative
				else
					dipole_derivatives = nil
				end
			else
				# When Hessian is read, read it into results
				Cuby::log.puts_debug("Reading hessian from file '#{@settings[:hessian_read]}'")
				f =  File.open(@settings[:hessian_read])
				@results.hessian = YAML.load(f)
				# Check if it is Matrix-like class (Hessian or Matrix, saved differently on difefrent rubies)
				unless @results.hessian.kind_of?(Matrix)
					Cuby::error ("Error reading hessian from file, it is not a Matrix-like object")
				end
				f.close
				# No dipole derivatives when Hessian is read from a file
				dipole_derivatives = nil
			end

			# Vibrational and thermodynamic analysis	
			energy = 0.0
			@vibrations = Vibrations.from_hessian(@geometry, @results.hessian, energy, @settings, dipole_derivatives)
			@thermodynamics = Thermodynamics.new(@settings, @geometry, @vibrations.frequencies, nil)

		end
		return @results
	end

	def transformation_normal_to_cart(vibrations)
		# Transformation matrix from dimensionless normal coordinates to cartesian displacement
		cols = []
		vibrations.normal_modes.each_index{|i|
			next if vibrations.frequencies[i] == 0.0
			factor = (CM2KCALMOL *  vibrations.frequencies[i]/vibrations.force_constants[i])**0.5
			vector = vibrations.normal_modes[i]
			cols << (vector * factor).to_a

		}
		m = Matrix.columns(cols)

		#--------------------------------------------------
		# # Test
		# gv = @geometry.to_vector
		# g = @geometry.deep_copy
		# v = Vector.zero(vibrations.normal_modes.size - 6)
		# v[0] = -1.0
		# g.update_from_vector!(gv + m*v)
		# g.write_xyz
		# v[0] = 0.0
		# g.update_from_vector!(gv + m*v)
		# g.write_xyz
		# v[0] = 1.0
		# g.update_from_vector!(gv + m*v)
		# g.write_xyz
		#-------------------------------------------------- 

		return m
	end

	def one_dimensional_scans(vibrations)
		Cuby::log.puts_v_only(:normal,"")
		Cuby::log.puts_v(:brief, "Generating scans in normal modes, force constant k in kcal/mol/A^2,\ndisplacement interval as multipliers of the normal mode:\n")
		Cuby::log.puts_v(:brief, "   mode   k                from                 to")
		vibrations.normal_modes.each_index{|i|
			next if vibrations.frequencies[i] == 0.0
			# Force constant
			k = vibrations.force_constants[i]
			case @settings[:freq_mode_scan_e_unit]
			when :fundamental_frequency
				# Max displacement: calculated so that the energy is x * v0
				# E = 1/2 * k * x^2
				e_max = @settings[:freq_mode_scan_e] * vibrations.frequencies[i] * CM2KCALMOL
			when :cm
				e_max = @settings[:freq_mode_scan_e] * CM2KCALMOL
			when :kcal
				e_max = @settings[:freq_mode_scan_e]
			end

			x_max = (2.0 * e_max.abs / k)**0.5

			Cuby::log.puts_v(:brief, sprintf("   %-5d%10.3f%20.8f%20.8f",i, k, -x_max, x_max))

			# Do the scan
			f = File.open("mode#{"%03d" % i}.xyz", "w+")
			step = 2.0 / (@settings[:freq_mode_scan_n] - 1)
			d = -1.0
			newg = @geometry.deep_copy
			gvec0 = @geometry.to_vector
			while d <= 1.0
				x = x_max * d
				gvec = gvec0 + vibrations.normal_modes[i] * x
				newg.update_from_vector!(gvec)	
				# q: value in dimensionless coordinates
				# In these coordinates, E[cm-1] = 0.5 * f0[cm-1]  * q**2
				q = (k * x**2 / CM2KCALMOL / vibrations.frequencies[i])**0.5 * x/x.abs
				newg.write_xyz(:file => f, :second_line => sprintf("scan: %.6f", q))
				d += step
			end
			f.close
		}
	end

	def multi_dimensional_scan(vibrations, modes)
		Cuby::log.puts_v_only(:normal,"")
		Cuby::log.puts_v(:brief, "Multi-dimensional scan in #{modes.size} modes")
		# Get displacements for each scanned mode
		displacements = ArrayOfArrays.new
		q_coords = ArrayOfArrays.new

		if @settings.set?(:freq_mode_scan_grids)
			modes.each_index{|mi|
				i = modes[mi]

				# Input is in dimensionless coordinates
				# In these coordinates, E[cm-1] = 0.5 * f0[cm-1]  * q**2
				# Calculate conversion factor to kcal/mol / A system
				conv_factor = (vibrations.frequencies[i] * CM2KCALMOL / vibrations.force_constants[i])**0.5

				item = @settings[:freq_mode_scan_grids][mi]
				disps = []
				q = []
				if item.class == String
					seq = Settings::Seq.from_string(item)
					seq.iterate{|d|
						disps << conv_factor * d
						q << d
					}
				elsif item.class == Array
					item.each{|d|
						disps << conv_factor * d
						q << d
					}
				else
					Cuby::error "Item of the freq_mode_scan_grids list must be either Array or a String"
				end

				displacements << disps
				q_coords << q
			}
		else
			# The same energy-based criteria for all modes
			modes.each{|i|
				k = vibrations.force_constants[i]
				# Force constant
				k = vibrations.force_constants[i]
				case @settings[:freq_mode_scan_e_unit]
				when :fundamental_frequency
					# Max displacement: calculated so that the energy is x * v0
					# E = 1/2 * k * x^2
					e_max = @settings[:freq_mode_scan_e] * vibrations.frequencies[i] * CM2KCALMOL
				when :cm
					e_max = @settings[:freq_mode_scan_e] * CM2KCALMOL
				when :kcal
					e_max = @settings[:freq_mode_scan_e]
				end
				x_max = (2.0 * e_max.abs / k)**0.5
				step = 2.0 / (@settings[:freq_mode_scan_n] - 1)
				d = -1.0
				disps = []
				q = []
				while d <= 1.0
					x = x_max * d
					disps << x
					q << (k / CM2KCALMOL / vibrations.frequencies[i])**0.5 * x
					d += step
				end
				displacements << disps
				q_coords << q
			}
		end
		total = q_coords.no_of_combinations
		Cuby::log.puts_v(:brief, "Number of points: #{total}")
	
		# Print the grid
		Cuby::log.puts_v_only(:normal,"")
		Cuby::log.puts_v_only(:normal,"Displacement grid (in dimensionless coordinates)")
		modes.each_index{|i|
			Cuby::log.puts_v_only(:normal,"mode #{modes[i]}: " +  q_coords[i].map{|x| "%.6f" % x}.join("\t"))
		}

		# Iterate over all combinations
		Cuby::log.puts_v_only(:normal,"")
		Cuby::log.puts_v_only(:normal,"Generating points:")
		m = modes.size
		newg = @geometry.deep_copy
		gvec0 = @geometry.to_vector
		if total < @settings[:freq_mode_scan_batchsize]
			f = File.open("multidimensional_scan.xyz", "w+")
			Cuby::log.puts_v_only(:normal,"  #{0} - #{total - 1} / #{total}")
		else
			dec = (total-1).to_s.size
			batch_s = 0
			batch_e = @settings[:freq_mode_scan_batchsize] - 1
			f = File.open(sprintf("multidimensional_scan_%0#{dec}d-%0#{dec}d.xyz", batch_s, batch_e), "w+")
			Cuby::log.puts_v_only(:normal,"  #{batch_s} - #{batch_e} / #{total}")
		end
		i = 0
		q_coords.iterate_all_combinations{|v|
			# v is the array of indices of displacement in each mode
			vec= gvec0.deep_copy
			s = ""
			m.times{|mi|
				ii = modes[mi] # index of the mode in the list
				vec += vibrations.normal_modes[ii] * displacements[mi][v[mi]]
				s += sprintf(" %.6f",q_coords[mi][v[mi]])
			}
			newg.update_from_vector!(vec)
			newg.write_xyz(:file => f, :second_line => "scan: "+s)

			# New batch
			if ((i+1) % @settings[:freq_mode_scan_batchsize] == 0 ) && (i+1 != total)
				batch_s += @settings[:freq_mode_scan_batchsize]
				batch_e += @settings[:freq_mode_scan_batchsize]
				batch_e = total - 1 if batch_e >= total
				f.close
				f = File.open(sprintf("multidimensional_scan_%0#{dec}d-%0#{dec}d.xyz", batch_s, batch_e), "w+")
				Cuby::log.puts_v_only(:normal,"  #{batch_s} - #{batch_e} / #{total}")
			end
			i += 1
			
		}
		f.close
	end
end
