
class Hybridization < Array
	attr_reader :gradient # Gradient of the continuous hybridiation

	def initialize(geometry)
		@geometry = geometry
		@gradient = nil
	end

	#=======================================================================
	# Integer hybridization indices from connectivity
	#=======================================================================

	def evaluate_from_connectivity(connectivity = nil)
		#: Integer hybridization, bond counts from connectivity

		unless connectivity
			# Calculate connectivity if not provided
			connectivity = Connectivity.new(@geometry)
		end
		@geometry.each_index {|i|
			self[i] = connectivity.hybridization(i)
		}
	end

	#=======================================================================
	# Grimme's continuous hybridization
	#=======================================================================

	# Constants
	GRIMME_NONMETALS  = [:H, :He, :C, :N, :O, :F, :Ne, :P, :S, :Cl, :Ar, :Br, :Kr, :I, :Xe, :Hg, :Rn]
	### Possible error in Grimme's code - Hg listed as non-metal
	GRIMME_HALFMETALS = [:As, :Se, :At]

	def precalculate_grimme_radii(covalent_radii, scale_radii)
		radii = {}
		@geometry.elements_in_system.each{|element|
			# Get covalent radius
			if covalent_radii.kind_of?(Hash)
				# Hash - table of radii
				r = covalent_radii[element]
			elsif covalent_radii.kind_of?(Symbol)
				# Symbol - name of radii set available in PeriodicTable
				r = PeriodicTable.covalent_radius(element, covalent_radii)
			else
				raise "Covalent radii has to be provided as a list or name of set in PeriodicTable"
			end
			# Scaling according to Grimme's paper
			if scale_radii
				# Scale for metals
				if GRIMME_NONMETALS.include?(element)
					r *= 1.0
				elsif GRIMME_HALFMETALS.include?(element)
					r *= 0.95
				else
					r *= 0.9
				end
			end
			# Save
			radii[element] = r
		}
		return radii
	end

	def grimme_atom(i, radii,  grad = false)
		atom_i = @geometry.at(i)


		hyb = 0.0
		if grad
			gradient = Coordinate.new
		end
		@geometry.each_index {|j|
			next if i == j
			atom_j = @geometry.at(j)

			r = atom_i.distance(atom_j)
			r_thr = radii[atom_i.element] + radii[atom_j.element]

			cn = 1.0 / (1.0 + Math::exp(-16.0 * (4.0/3.0 * r_thr / r - 1.0)))
			hyb += cn

			if grad
				deriv = -(64.0*(r_j+r_i) * Math::exp(64.0*r_thr/(3.0*r) + 16.0)) / 
					(3.0*r**2*(Math::exp(64*r_thr/(3*r)) + Math::exp(16))**2)
				gradient += deriv * (atom_i - atom_j) / r
			end
		}

		if grad
			return [hyb, gradient]
		else
			return hyb
		end
	end

	def evaluate_grimme(covalent_radii = :PA2009, scale = true, grad = false)
		#: From Grimme's DFT-D3 paper
		### Add full reference when available

		# Precalculate radii, once per geometry is enough
		unless @radii
			@radii = precalculate_grimme_radii(covalent_radii, scale)
		end

		# Fill the array(s)
		if grad
			@gradient = []
			@geometry.each_index {|i|
				self[i], @gradient[i] = grimme_atom(i, @radii, grad)
			}
		else
			@geometry.each_index {|i|
				self[i] = grimme_atom(i, @radii, grad)
			}
		end

		return nil
	end

	#=======================================================================
	# Grimme-like with finite cutoff
	#=======================================================================

	def grimme_mod_atom(i, radii, grad = false)
		atom_i = @geometry.at(i)

		hyb = 0.0
		if grad
			gradient = Coordinate.new
		end
		@geometry.each_index {|j|
			next if i == j
			atom_j = @geometry.at(j)
			r = atom_i.distance(atom_j)

			r_cut = radii[atom_i.element] + radii[atom_j.element]
			r1 = r_cut * 1.75

			# Above threshold - can be skipped
			next if r >= r1

			# Other options
			r0 = r_cut * 0.95
			if r <= r0
				cn = 1.0
			else
				x = (r - r0) /(r1 - r0)
				cn = 1.0 - (-20.0*x**7 + 70.0*x**6 - 84.0*x**5 + 35.0*x**4)
			end		
			hyb += cn

			if grad
				raise "No grad yet"
			end
		}

		if grad
			return [hyb, gradient]
		else
			return hyb
		end
	end

	def evaluate_grimme_mod(covalent_radii = :PA2009, scale = true, grad = false)
		# Modified version of Grimme's approach
		# Sigmoidal function with finite cutoff radius is used, not all pairs have
		# to be evaluated

		# Precalculate radii, once per geometry is enough
		unless @radii
			@radii = precalculate_grimme_radii(covalent_radii, scale)
		end

		# Fill the array(s)
		if grad
			@gradient = []
			@geometry.each_index {|i|
				self[i], @gradient[i] = grimme_mod_atom(i, @radii, grad)
			}
		else
			@geometry.each_index {|i|
				self[i] = grimme_mod_atom(i, @radii, grad)
			}
		end

		return nil
	end
end
