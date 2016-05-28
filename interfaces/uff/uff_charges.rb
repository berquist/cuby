class UffCharges
	
	C = 14.4
	
	#~ def coulomb(x)   # Coulomb integral for atom centres distance x; charge[e], distance[\circ A]
		#~ return C*(1-Math.exp(-x-0.5*x**2))/x if x != 0
		#~ return C
		#~ # return 1/(x+2.5)
	def coulomb(x, j)   # Coulomb integral for atom centres distance x; charge[e], distance[\circ A]
		return C*(1-Math.exp(-j/C*x-j**2*x**2/(2*C**2)))/x if x != 0
		return j
	end
	
	def initialize(geometry, charge)
		evaluate_charges(geometry, charge)
		Cuby.log.puts_debug "UFF charges generated:"
		geometry.each_index{|i|
			Cuby.log.puts_debug sprintf("  %3s   %8.3f", geometry.at(i).element.to_s, @charges[i])
		}
	end

	def [] (index)
		return @charges[index]
	end

	#=======================================================================
	# Atomic parameters
	#=======================================================================

	# Electronegativity, eV
	X = {
		:H => 4.528,
		:C => 5.343,
		:N => 6.899,
		:O => 8.741,
		:Na => 2.843,
		:Cl => 8.564,
	}

	# Idempotential, eV
	J = {
		:H => 13.8904,
		:C => 10.126,
		:N => 11.760,
		:O => 13.364,
		:Na => 4.592,
		:Cl => 9.892,
	}
	
	#=======================================================================
	# Private methods
	#=======================================================================

	def evaluate_charges(geometry, charge)
		@charges = []
		n = geometry.size
		b = []
		geometry.each_index{|i|
			@charges[i] = 0.0
			b[i] = -X[geometry[i].element]   # fills vector with electronegativities
		}
		b[n] = charge
		
		# possibly sparse matrix for huge systems
		
		a = Matrix.zero(n+1, n+1)
		n.times{|i|
			a[n,i] = 1
			a[i,n] = 1
			n.times{|j|
				a[i,j] = coulomb(((geometry[i].x-geometry[j].x)**2 + (geometry[i].y-geometry[j].y)**2 + (geometry[i].z-geometry[j].z)**2)**0.5, [J[geometry[i].element] , J[geometry[j].element]].min)   # porad spatne !!!, navic chybi omezeni naboju !!!
				puts coulomb(((geometry[i].x-geometry[j].x)**2 + (geometry[i].y-geometry[j].y)**2 + (geometry[i].z-geometry[j].z)**2)**0.5, (J[geometry[i].element] + J[geometry[j].element])/2)   if i != j
			}
		}
		
		puts "Charges:"
		q = (a.inverse * Vector.from_array(b)).to_a
		q.slice!(-1)
		puts q
		
		return nil
	end
	
end
