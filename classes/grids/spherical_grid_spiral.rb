class SphericalGrid
	attr_reader :coords

	#=======================================================================
	# Initialization
	#=======================================================================
	
	def initialize(r, n)
		@r = r
		@n = n
		@coords = nil
	end

	def build_coords(center = Coordinate[0,0,0])
		@coords = []
		iterate_xyz(center){|crd|
			@coord << crd
		}
		return true
	end

	#=======================================================================
	# Sphere properties
	#=======================================================================
	
	def surface
		return 4.0 * Math::PI * @r**2
	end
	
	def volume
		return 4.0/3.0 * Math::PI * @r**3
	end

	#=======================================================================
	# Iterators
	#=======================================================================
	
	def iterate_xyz(center = Coordinate[0,0,0])
		iterate_spherical{|theta, phi|
			yield SphericalGrid.spherical_to_cartesian(theta, phi, @r) + center
		}
		return @n
	end

	#=======================================================================
	# Helpers
	#=======================================================================
	
	def SphericalGrid.spherical_to_cartesian(theta, phi, r)
		x = r * Math.sin(theta) * Math.cos(phi)
		y = r * Math.sin(theta) * Math.sin(phi)
		z = r * Math.cos(theta)
		return Coordinate[x,y,z]
	end

	def SphericalGrid.cartesian_to_spherical(x,y,z)
		r = (x**2 + y**2 + z**2)**0.5
		theta = Math.acos(z/r)
		phi = Math.atan2(y,x)
		return theta, phi, r
	end
	
end

class SphericalGridSpiral < SphericalGrid

	def iterate_spherical
		# Saff and Kuijlaars algorithm
		# http://sitemason.vanderbilt.edu/page/hmbADS
		p = 0.5
		a = 1.0 - 2.0 * p / (@n - 3)
		b = p * (1.0 + @n) / (@n - 3)
	
		rk = 0.0

		# First point
		theta = Math::PI
		phi = 0.0
		yield(theta, phi)

		# All points in between
		(@n-2).times{|kk|
			k = kk + 2
			kb = a*k + b

			hk = -1.0 + 2.0 * (kb-1.0)/(@n-1.0)
			rk_prev = rk
			rk = (1.0 - hk**2)**0.5

			theta = Math.acos(hk)
			phi = (phi + 3.6 / @n**0.5 * 2.0 / (rk_prev + rk)) % 
			      (2.0 * Math::PI)
			yield(theta, phi)
		}

		# Last point
		theta = 0.0
		phi = 0.0
		yield(theta, phi)

		return nil
	end

end
