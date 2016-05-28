#===============================================================================
# Calculation results in a format read by gaussian
# Documentation: http://gaussian.com/g_tech/g_ur/k_external.htm
#===============================================================================

module GaussianExternalOutput

	def self.write(filename, results, what)
		f = File.open(filename, "w+")

		if what.include?(:energy)
			f.printf("%20.12e%20.12e%20.12e%20.12e\n", results.energy * KCAL2HARTREE, 0.0, 0.0, 0.0)
		end

		if what.include?(:gradient)
			results.gradient.each{|g|
				g = g * KCAL2HARTREE / ANGSTROM2BOHR
				f.printf("%20.12e%20.12e%20.12e\n", g.x, g.y, g.z)
			}
		end

		if what.include?(:hessian)
			# polarizability	    	Polar(I), I=1,6	    	3D20.12
			f.printf("%20.12e%20.12e%20.12e\n", 0.0, 0.0, 0.0)
			f.printf("%20.12e%20.12e%20.12e\n", 0.0, 0.0, 0.0)

			# dipole derivatives	    	DDip(I), I=1,9*NAtoms	    	3D20.12
			results.hessian.n.times{
				f.printf("%20.12e%20.12e%20.12e\n", 0.0, 0.0, 0.0)
			}

			# force constants	    	FFX(I), I=1,(3*NAtoms*(3*NAtoms+1))/2	    	3D20.12
			n = results.hessian.n
			count = 0
			n.times{|i|
				(i+1).times{|j|
					f.printf("%20.12e", results.hessian[i,j] * KCAL2HARTREE / ANGSTROM2BOHR / ANGSTROM2BOHR)
					count += 1
					f.puts if count % 3 == 0
				}
			}
			
		end

		f.close
	end

end
