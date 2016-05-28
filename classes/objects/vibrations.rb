################################################################################
#
# Class Vibrations
#
# Author: Jan Rezac
# Date created: 2009-08-24
# License: Cuby license
# Description: %short_description
# Status: Planned
#
################################################################################

#: More detailed description

class VibrationsData
	attr_accessor :energy			# Float
	attr_accessor :frequencies 		# Array of Floats
	attr_accessor :normal_modes		# Array of Vectors, in cartesian coordinates
	attr_accessor :normal_modes_mw		# Array of Vectors, in mass-weighted coordinates
	attr_accessor :reduced_mass		# Array of Floats
	attr_accessor :intensities		# Array of Floats
	attr_accessor :dipole_derivatives	#
	attr_accessor :force_constants		# Array of Floats

	#=======================================================================
	# Save/Load
	#=======================================================================

	def VibrationsData.load_yaml(file)
		if file.class == String
			file = File.open(file,"r")
			closefile = true
		end

		new = YAML.load(file)
		file.close if closefile
		return new
	end

	def save_yaml(file)
		if file.class == String
			file = File.open(file,"w+")
			closefile = true
		end

		file.puts(self.to_yaml)
		file.close if closefile
	end
end

class Vibrations

	#=======================================================================
	# Initialize
	#=======================================================================

	def initialize(frequencies_data)
		#: for internal use only
		@data = frequencies_data
	end

	def Vibrations.from_file(file)
		return Vibrations.new(VibrationsData.load_yaml(file))
	end

	def Vibrations.from_hessian(geometry, hessian, energy, settings, dipole_derivatives)
		# Change atom masses to pure isotopes, but keep the user-defined values
		geometry.each{|atom|
			atom.properties[:mass] = PeriodicTable.common_isotope(atom.element).mass unless atom.properties[:mass]
		}

		#if true
		if false
			return Vibrations.from_hessian_new(geometry, hessian, energy, settings, dipole_derivatives)
		else
			return Vibrations.from_hessian_old(geometry, hessian, energy, settings, dipole_derivatives)
		end
	end

	def Vibrations.from_hessian_new(geometry, hessian, energy, settings, dipole_derivatives)

		# Translation / rotation removal
		projector = Hessian.transrot_removal_projector(geometry)
		hessian = projector.transpose * hessian * projector

		# Mass weighting
		mass_matrix = Hessian.mass_unweight_matrix(geometry)
		mass_matrix_inv = Hessian.mass_weight_matrix(geometry)
		hessian_mw = mass_matrix_inv * hessian * mass_matrix_inv

		# Diagonalization
		eigenvectors, real, imag = hessian_mw.eigensystem # eigenvectors == L matrix
		Matrix.sort_eigensystem!(eigenvectors, real, imag, :descending_real)


		# Eigenvalues converted to ferquencies
		# Eigenvalues are in 1/ps^2
		frequencies = []
		puts "eigenvalues"
		real.each{|x| 
			frequencies << x**0.5 / PS2S / 2.0 / Math::PI *  HZ2CM
			puts x**0.5 / PS2S / 2.0 / Math::PI *  HZ2CM
			#puts (x / PS2S**2 / 4.0 / Math::PI**2 / SPEED_OF_LIGHT_SI**2)**0.5 / 100
		}
		puts

		# Normal modes in cartesian coordinates
		modes_in_cart = (eigenvectors.transpose * mass_matrix).inverse # this keeps the translation/rotation modes
		# mass_matrix_inv * projector * eigenvectors # this is equivalent for vibrations, but zeroes out translation and rotation normal vectors

		puts modes_in_cart
		puts
		puts mass_matrix_inv * eigenvectors
		modes_in_cart = mass_matrix_inv * eigenvectors

		# reduced masses
		puts "reduced masses"
		redmass = []
		modes_in_cart.n.times{|i|
			redmass[i] = modes_in_cart.column_as_vector(i).abs**-2 * UNIT2GMOL
		}
		puts redmass
		puts

		# Normalize the cartesian normal vectors
		cart_vec = []
		modes_in_cart.n.times{|i|
			# using calculated normalization constants
			#cart_vec[i] = modes_in_cart.column_as_vector(i) * (redmass[i]/UNIT2GMOL)**0.5
			# simpler equivalent
			cart_vec[i] = modes_in_cart.column_as_vector(i).normalize
		}

		puts "normal mode0"
		puts cart_vec[0]


		# Force constants in normal modes
		fcmat = eigenvectors.transpose * hessian_mw * eigenvectors
		fcarray = []
		fcmat.each_diagonal{|x| fcarray << x}
		fcvec = Vector.from_array(fcarray)

		puts "force constants"
		#puts (fcvec).to_a
		# to get the FC, it must be divided by the size of the normal mode ^2 (it is the reduced mass)
		puts fcvec[0] / modes_in_cart.column_as_vector(0).abs**2

		# Another way of getting the FCs
		puts (redmass[0]/UNIT2GMOL) * (2.0 * Math::PI * frequencies[0] * CM2HZ * PS2S)**2

		# Another?
		fc_mwc = real[0,0] * (redmass[0]/UNIT2GMOL)
		puts fc_mwc

		# Do 1kcal/mol E displacement
		e_max = 1.23456
		k = fc_mwc
		x_max = (2.0 * e_max / k)**0.5
		puts "x_max #{x_max}"

		puts cart_vec[0]
		vec = cart_vec[0]
		puts vec.abs
		vec = vec * x_max
	
		f = File.open("test.xyz", "w+")
		geometry.write_xyz(:file => f)
		gvec0 = geometry.to_vector
		geometry.update_from_vector!(gvec0 + vec)
		geometry.write_xyz(:file => f)
		f.close

		exit
	end

	def Vibrations.from_hessian_old(geometry, hessian, energy, settings, dipole_derivatives)
		# translation / rotation removal
		projector = Hessian.transrot_removal_projector(geometry)
		hessian = projector.transpose * hessian * projector

		# mass weight hessian
		weighting_matrix = Hessian.mass_weight_matrix(geometry)
		hessian = weighting_matrix * hessian * weighting_matrix

		# Diagonalization
		vectors, real, imag = hessian.eigensystem

		# Sort eigensystem
		Matrix.sort_eigensystem!(vectors, real, imag, :descending_real)

		# Convert eigenvalues to arrays
		real = real.to_a.flatten
		imag = imag.to_a.flatten


		# Process the normal vectors, saving following data:
		vectors_mw = []	# Mass-weighted normal vectors
		vectors_c = []	# Cartesian vectors
		redmass = [] 	# Reduced masses, in g/mol

		vectors.n.times{|v|
			# Save columns as mass-weighted normal vectors
			vectors_mw = vectors.column_as_vector(v)
			# Get normal vectors in cartesian coordinates (undo mass weighting), get reduced masses
			vectors.m.times{|i|
				vectors[i,v] = vectors[i,v] / (geometry[i/3].mass**0.5)
			}
			redmass << 1.0/vectors.column_as_vector(v).abs**2 * UNIT2GMOL
			# Normalize normal modes
			f = vectors.column_as_vector(v).abs
			vectors.m.times{|i|
				vectors[i,v] = vectors[i,v] / f
			}
			# Save the cartesian vector
			vectors_c[v] = vectors.column_as_vector(v)
		}


		# Intensities
		if dipole_derivatives
			# Three derivative vectors
			ddx = Vector.of_size(dipole_derivatives.size)
			ddy = Vector.of_size(dipole_derivatives.size)
			ddz = Vector.of_size(dipole_derivatives.size)
			dd = [ddx, ddy, ddz]
			dipole_derivatives.each_index{|i|
				3.times{|j|
					dd[j][i] = dipole_derivatives[i][j]
				}
			}

			# apply trans/rot removal to dipole derivatives vector
			3.times{|i|
				dd[i] = (projector * dd[i]).to_vector
			}

			# Transform to normal coordinates
			# Dot product of matrix column with the derivative
			3.times{|i|
				dd[i] = (vectors.transpose * dd[i]).to_vector
			}
			# Get squares of absolute values
			dipole_derivatives.each_index{|i|
				dipole_derivatives[i] = (dd[0][i]**2 + dd[1][i]**2 + dd[2][i]**2)
			}

			# Save as intensities
			intensities = dipole_derivatives.to_a
		else
			intensities = Array.new(real.size, 1.0)
		end

		# Calculate frequencies
		freqs = []
		real.each_index {|i|
			l = Complex(real[i],imag[i])
			l = 0.0 if l.abs < 1e-12
			sign = l / l.abs if l != 0.0
			l2 = l.abs
			if l == 0
				f = 0
			else
				f = l2**0.5 / 2 / Math::PI * sign
				f = f.real
				f /= 1e-12 # and convert it from 1/ps to 1/s
				f *= HZ2CM
			end
			f = 0.0 if f.abs < 1.0e-4
			freqs[i] = f
			if f == 0.0
				intensities[i] = 0.0
			end
		}

		# Calculate and convert force constants
		force_constants = []
		freqs.each_index{|i|	
			unless freqs[i] == 0.0
				# Force constant
				k = 4.0 * Math::PI**2 * redmass[i] * GMOL2UNIT * MASSUNIT2KG * (freqs[i] * CM2HZ)**2 / 100 # miliDyne/A
				# Unit conversion to Cuby units
				k = k /1000 # Dyne/A
				k = k * 1.0e-5 # N/A
				k = k * NEWTON2AU * AU2CUBYFORCE # kcal/mol/A^2
			else
				k = nil
			end
			# Save
			force_constants[i] = k
		}

		# Build the data structure
		data = VibrationsData.new
		data.energy = energy
		data.frequencies = freqs
		data.normal_modes = vectors_c
		data.normal_modes_mw = vectors_mw
		data.reduced_mass = redmass
		data.intensities = intensities
		data.dipole_derivatives = dipole_derivatives if dipole_derivatives
		data.force_constants = force_constants
		
		return Vibrations.new(data)
	end

	#--------------------------------------------------
	# def Vibrations.from_hessian2(geometry, hessian, energy, settings, gradient, dipole_derivatives)
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Partial hessian
	# 	list_negated = geometry.atomlist_from_selection("%not(" + settings[:partial_hessian] + ")")
	# 	list_negated.each{|at_i|
	# 		geometry[at_i].mass = 1.0e8
	# 	}
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# translation / rotation removal
	# 	#projector = Hessian.transrot_removal_projector_nomw(geometry)
	# 	#hessian = projector.transpose * hessian * projector
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# build hessian for the region
	# 	list = geometry.atomlist_from_selection(settings[:partial_hessian])
	# 	region = geometry.geometry_from_list(list)
	# 	r_hessian = Matrix.zero(list.size * 3)
	# 	r_grad = Gradient.new
	# 	r_dip = [] if dipole_derivatives
	# 	list.each_index{|i|
	# 		r_grad << gradient[list[i]]
	# 		if dipole_derivatives
	# 			r_dip << dipole_derivatives[list[i]*3]
	# 			r_dip << dipole_derivatives[list[i]*3+1]
	# 			r_dip << dipole_derivatives[list[i]*3+2]
	# 		end
	# 		list.each_index{|j|
	# 			3.times{|a|
	# 				3.times{|b|
	# 					r_hessian[3*i+a,3*j+b] = hessian[3*list[i]+a,3*list[j]+b]
	# 				}
	# 			}
	# 		}
	# 	}
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Trans/rot removal - Hessian
	# 	projector = Hessian.transrot_removal_projector_nomw(region)
	# 	r_hessian = projector.transpose * r_hessian * projector
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Reconstruct full hessian
	# 	hessian = Matrix.diagonal(geometry.size * 3, 1.0e-8)
	# 	#hessian = Matrix.zero(geometry.size * 3)
	# 	(list.size * 3).times{|i|
	# 	(list.size * 3).times{|j|
	# 		fi = list[i/3]*3 + i % 3
	# 		fj = list[j/3]*3 + j % 3
	# 		hessian[fi,fj] = r_hessian[i,j]
	# 	}
	# 	}
	#-------------------------------------------------- 


	#--------------------------------------------------
	# 	# trans/rot removal: gradient
	# 	g = r_grad.to_vector
	# 	g = (projector * g).to_vector
	# 	# Gradient control
	# 	if g.max_abs > 0.6
	# 		puts
	# 		puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	# 		puts "Max. gradient: #{sprintf('%5.2f',maxgrad)} indicates that the structure is not a stationary point"
	# 		puts "Calculated frequencies have no meaning"
	# 		puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
	# 		puts
	# 	end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Trans/rot removal - dipole derivatives
	# 	if dipole_derivatives
	# 		# Build vectors
	# 		ddx = Vector.of_size(r_dip.size)
	# 		ddy = Vector.of_size(r_dip.size)
	# 		ddz = Vector.of_size(r_dip.size)
	# 		r_grad.size.times{|i|
	# 			ddx[i] = r_dip[i].x
	# 			ddy[i] = r_dip[i].y
	# 			ddz[i] = r_dip[i].z
	# 		}
	# 		ddx = (projector * ddx).to_vector
	# 		ddy = (projector * ddy).to_vector
	# 		ddz = (projector * ddz).to_vector
	# 	end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Reconstruct dipole derivatives
	# 	if dipole_derivatives
	# 		ddxf = Vector.zero(geometry.size * 3)
	# 		ddyf = Vector.zero(geometry.size * 3)
	# 		ddzf = Vector.zero(geometry.size * 3)
	# 		(list.size * 3).times{|i|
	# 			fi = list[i/3]*3 + i % 3
	# 			ddxf[fi] = ddx[i]
	# 			ddyf[fi] = ddy[i]
	# 			ddzf[fi] = ddz[i]
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 		}
	# 		ddx = ddxf
	# 		ddy = ddyf
	# 		ddz = ddzf
	# 	end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# mass weight hessian
	# 	weighting_matrix = Hessian.mass_weight_matrix(geometry)
	# 	hessian = weighting_matrix * hessian * weighting_matrix
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Diagonalization
	# 	vectors, real, imag = hessian.eigensystem
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Sort eigensystem
	# 	Matrix.sort_eigensystem!(vectors, real, imag, :descending_real)
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Get normal vectors in cartesian coordinates (undo mass weighting), get reduced masses
	# 	redmass = []
	# 	vectors.n.times{|v|
	# 		vectors.m.times{|i|
	# 			#vectors[i,v] = vectors[i,v] / (geometry[i/3].mass**0.5)
	# 		}
	# 		redmass << 1.0/vectors.column_as_vector(v).abs**2 * UNIT2GMOL
	# 	}
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Normalize normal modes
	# 	vectors.n.times{|v|
	# 		f = vectors.column_as_vector(v).abs
	# 		vectors.m.times{|i|
	# 			vectors[i,v] = vectors[i,v] / f
	# 		}
	# 	}
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Intensities
	# 	if dipole_derivatives
	# 		# Three derivative vectors
	# 		dd = [ddx, ddy, ddz]
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 		# Transform to normal coordinates
	# 		# Dot product of matrix column with the derivative
	# 		3.times{|i|
	# 			dd[i] = (vectors.transpose * dd[i]).to_vector
	# 		}
	# 		# Get squares of absolute values
	# 		dd[0].each_index{|i|
	# 			dipole_derivatives[i] = (dd[0][i]**2 + dd[1][i]**2 + dd[2][i]**2)
	# 		}
	# 	end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Convert everything to arrays
	# 	sorted_r = real.to_a.flatten
	# 	sorted_i = imag.to_a.flatten
	# 	vectors_c = []
	# 	vectors.n.times{|i|
	# 		vectors_c[i] = vectors.column_as_vector(i)
	# 	}
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Intensities
	# 	if dipole_derivatives
	# 		intensities = dipole_derivatives.to_a
	# 	else
	# 		intensities = Array.new(sorted_r.size, 1.0)
	# 	end
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Calculate frequencies
	# 	freqs = []
	# 	sorted_r.each_index {|i|
	# 		l = Complex(sorted_r[i],sorted_i[i])
	# 		l = 0.0 if l.abs < 1e-12
	# 		sign = l / l.abs if l != 0.0
	# 		l2 = l.abs
	# 		if l == 0
	# 			f = 0
	# 		else
	# 			f = l2**0.5 / 2 / Math::PI * sign
	# 			f = f.real
	# 			f /= 1e-12 # and convert it from 1/ps to 1/s
	# 			f *= HZ2CM
	# 		end
	# 		f = 0.0 if f.abs < 1.0e-4
	# 		freqs[i] = f
	# 		if f == 0.0
	# 			intensities[i] = 0.0
	# 		end
	# 	}
	#-------------------------------------------------- 

	#--------------------------------------------------
	# 	# Build the data structure
	# 	data = VibrationsData.new
	# 	data.energy = energy
	# 	data.frequencies = freqs
	# 	data.normal_modes = vectors_c
	# 	data.reduced_mass = redmass
	# 	data.intensities = intensities
	# 	
	# 	return Vibrations.new(data)
	# end
	#-------------------------------------------------- 

	#=======================================================================
	# Data access
	#=======================================================================
	
	def frequencies
		return @data.frequencies
	end

	def normal_modes
		return @data.normal_modes
	end

	def reduced_mass
		return @data.reduced_mass
	end

	def energy
		return @data.energy
	end

	def intensities
		return @data.intensities
	end

	def dipole_derivatives
		return @data.dipole_derivatives
	end

	def force_constants
		return @data.force_constants
	end
	
	#=======================================================================
	# Input/Output
	#=======================================================================

	def to_file(file)
		@data.save_yaml(file)
	end

	def write_molden_input(geometry, settings, filename = "molden.input")
		f = File.open(filename, "w+")
		f.puts "[Molden Format]"
		f.puts "[Title]"
		f.puts "title"

		if settings.set?(:partial_hessian)
			partial = true
			list = geometry.atomlist_from_selection(settings[:partial_hessian])
			partial_geo = !settings[:partial_hessian_full_geometry]
		end

		f.puts "[Atoms] Angs"
		if partial && partial_geo
			list.each {|i|
				count = 1
				a = geometry[i]
				f.printf("%4s %4d %4d %12.6f %12.6f %12.6f\n",a.element.to_s.downcase, count, a.proton_number, a.x, a.y, a.z)
				count += 1
			}
		else
			geometry.each_index {|i|
				a = geometry[i]
				f.printf("%4s %4d %4d %12.6f %12.6f %12.6f\n",a.element.to_s.downcase, i+1, a.proton_number, a.x, a.y, a.z)
			}
		end

		f.puts "[FREQ]"
		frequencies.each_index{|ii| 
			i = frequencies.size - 1 - ii
			if partial && frequencies[i] != 0.0 || !partial
				f.printf("%12.6f\n",frequencies[i])
			end
		}

		f.puts "[FR-COORD]"
		if partial && partial_geo
			list.each {|i|
				a = geometry[i]
				f.printf("%4s %12.6f %12.6f %12.6f\n",a.element.to_s.downcase, a.x * ANGSTROM2BOHR, a.y * ANGSTROM2BOHR, a.z * ANGSTROM2BOHR)
			}
		else
			geometry.each_index {|i|
				a = geometry[i]
				f.printf("%4s %12.6f %12.6f %12.6f\n",a.element.to_s.downcase, a.x * ANGSTROM2BOHR, a.y * ANGSTROM2BOHR, a.z * ANGSTROM2BOHR)
			}
		end

		f.puts "[FR-NORM-COORD]"
		count = 1
		frequencies.each_index{|ii|
			i = frequencies.size - 1 - ii

			if partial && frequencies[i] != 0.0 || !partial
				f.puts " vibration #{count}"
				count += 1
				if partial && partial_geo
					list.size.times {|j|
						3.times {|k|
							f.printf(" %12.6f",normal_modes[i][list[j]*3+k])
						}
						f.puts
					}
				else
					geometry.size.times {|j|
						3.times {|k|
							f.printf(" %12.6f",normal_modes[i][j*3+k])
						}
						f.puts
					}
				end
			end
		}

		f.puts "[INT]"
		frequencies.each_index{|ii| 
			i = frequencies.size - 1 - ii
			if partial && frequencies[i] != 0.0 || !partial
				f.printf("%12.6f\n",intensities[i])
			end
		}

		f.close
	end

	#=======================================================================
	# Calculations
	#=======================================================================
	
	def zpve
		zpve = 0.0
		n = 6
		n = 5 if frequencies.size == 6
		(frequencies.size - n).times {|i|
			zpve += frequencies[i] * 0.5 * CM2KCALMOL
		}
		return zpve
	end

	#=======================================================================
	# Printing
	#=======================================================================

	def print(settings = {})
		puts "Frequencies (cm^-1):"
		frequencies.each_index {|i|
			printf("%-6d f= %10.3f",i, frequencies[i])
			printf("        intensity = %10.3f", intensities[i]) if dipole_derivatives
			#printf("        red. mass = %10.3f", reduced_mass[i]) if @settings[:print] =~ /redmass/i
			puts
		}

		# Print ZPVE for each mode
		if settings[:freq_print].include?(:mode_zpve)
			puts
			puts "ZPVE per mode (kcal/mol)"
			frequencies.each_index {|i|
				printf("%-6d f= %10.6f\n",i, frequencies[i] * 0.5 * CM2KCALMOL)
			}
		end

		# Print reduced mass of the modes
		if settings[:freq_print].include?(:reduced_mass)
			puts
			puts "Reduced masses of the modes (g/mol)"
			frequencies.each_index {|i|
				printf("%-6d f= %10.6f\n",i, reduced_mass[i])
			}
		end

		# Print force constants
		if settings[:freq_print].include?(:force_constants)
			puts
			puts "Force constants of the modes (kcal/mol/A^2)"
			frequencies.each_index {|i|
				if  force_constants[i]
					printf("%-6d f= %10.6f\n",i, force_constants[i])
				end
			}
		end

		# Print ZPVE
		puts
		puts "ZPVE = #{"%.6f"%zpve} kcal/mol"

		write_spectrum("spectrum_ir.txt")
	end

	def write_spectrum(file)
		if file.class == String
			file = File.open(file, "w+")
			closefile = true
		else
			closefile = false
		end

		if true
			file.puts "0.0\t0.0"
			frequencies.each_index {|i|
				intensity = 1.0
				intensity = intensities[i] if dipole_derivatives
				file.puts "#{frequencies[i]}\t0.0"
				file.puts "#{frequencies[i]}\t#{intensity}"
				file.puts "#{frequencies[i]}\t0.0"
		}
		end

		file.close if closefile
	end

end
