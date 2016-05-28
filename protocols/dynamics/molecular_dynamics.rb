require "protocols/dynamics/thermostat.rb"

class MolecularDynamics
	attr_accessor :cycle
	attr_accessor :energy
	attr_accessor :gradient
	attr_reader   :dt, :dt2 # Timestep, square is stored separately for efficiency

	# Geometry and selection of atoms (from the full geoemtry) MD is applied to
	attr_reader   :geometry_full
	attr_reader   :geometry
	attr_reader   :atomlist

	# Handles
	attr_accessor :print_step_proc

	def initialize(geometry, settings)
		# Keep settings locally
		@settings = settings

		# Get geometry on which MD is performed
		if @settings[:md_region] == "%all()"
			# The whole system is simulated
			@geometry = geometry
		else
			# Only a part of the system is moving
			@atomlist = geometry.atomlist_from_selection(@settings[:md_region])
			@geometry = geometry.geometry_from_list(@atomlist)
			@geometry_full = geometry

			# Switch off trans/rot removal
			if @settings[:remove_translation]
				@settings[:remove_translation] = false
				Cuby::warning "Freezer active, switching off translation removal"
			end
			if @settings[:remove_rotation]
				@settings[:remove_rotation] = false
				Cuby::warning "Freezer active, switching off rotation removal"
			end
		end

		# Random seed
		srand(@settings[:random_seed]) if @settings.set?(:random_seed)

		# Timestep
		self.dt = @settings[:timestep] # in ps

		# Initialize velocities
		initialize_velocities

		# Generate initial velocities
		case @settings[:velocities]
		when :zero
			# Nothing to be done
		when :random
			generate_maxwell_velocities(@settings[:init_temp])
		when :read
			if @settings[:velocities_file] =~ /.xyz[v]?$/
				g2 = Geometry.new
				g2.read_xyz(:file => @settings[:velocities_file], :velocities => true)
				if @settings[:md_region] != "%all()"
					gtemp = g2.geometry_from_list(@atomlist)
					g2 = gtemp
				end
				@geometry.each_index{|i|
					set_atom_velocity(i, g2[i].properties[:velocity])
				}
			else
				Cuby::error("Currently, velocities can be read only from .xyz file")
			end
		end

		# Set up thermostat
		case @settings[:thermostat]
		when :none
			extend Thermostat::None
			Cuby::error ("Keyword temperature_target can not be used without thermostat") if settings.set?(:temperature_target)
		when :berendsen
			extend Thermostat::Berendsen
		when :andersen
			extend Thermostat::Andersen
		when :"nose-hoover"
			extend Thermostat::NoseHoover
		when :bussi
			extend Thermostat::Bussi
		end
		thermostat_setup
	end

	def run
		# Initialize cycle counter
		@cycle = 0
		maxcycles = @settings[:maxcycles]
		while @cycle < maxcycles
			# Calculate
			@energy, @gradient = yield()

			# Call the print proc if it exists
			@print_step_proc.call unless @print_step_proc == nil

			# Integrate
			integrator

			# Fix COM translation/rotation
			remove_com_translation if @settings[:remove_translation]
			remove_com_rotation if @settings[:remove_rotation]
	
			# Interpolate and set thermostat temperature
			if @settings.set?(:temperature_target)
				thermostat_set_t(interpolate_temperature(@cycle))
			end

			# Thermostat
			e_kin, t = kinetic_and_temperature
			thermostat(t)

			@cycle += 1
		end
		@cycle -= 1
	end

	#=======================================================================
	# Access to attributes
	#=======================================================================
	
	def dt=(value)
		#: Sets the timestep, updating attribute dt2 (dt^2) as well.
		@dt2 = value**2
		@dt = value
	end

	#=======================================================================
	# Velocities
	#=======================================================================
	# The Verlet algorithm does not use velocities, but stores coordinates
	# from the previous step (variable @oldcoord)

	def initialize_velocities
		#: Sets atom velocities to zero
		@oldcoord = []
		@geometry.each_index {|i|
			# Copy coordinates of the atom to the t-1 coordinate
			@oldcoord[i] = Coordinate.from_coordinate(@geometry[i])
		}
	end

	def atom_velocity(i)
		return (@geometry[i] - @oldcoord[i]) / @dt
	end

	def set_atom_velocity(i, velocity)
		@oldcoord[i] = @geometry[i] - velocity * @dt
	end

	def velocities
		velo = []
		@geometry.each_index{|i|
			velo << atom_velocity(i)
		}
		return velo
	end

	def velocities_as_vector
		v = Vector.zero(@geometry.size * 3)
		@geometry.each_index{|i|
			crd = atom_velocity(i)
			3.times{|c|
				v[i*3+c] = crd[c]
			}
		}
		return v
	end

	def set_velocities_from_vector!(vector)
		@geometry.each_index{|i|
			atom_velo = Coordinate.new
			3.times{|c|
				atom_velo[c] = vector[i*3+c]
			}
			set_atom_velocity(i, atom_velo)
		}
		return nil
	end

	#=======================================================================
	# Generation of random velocities
	#=======================================================================

	def rand_gauss # => float
		#: Generator of random numbers with gaussian distribution. Currently,
		#: we use quick n' dirty algorithm, better one can be found in
		#: module mod_extend_math.
		return rand + rand + rand + rand + rand + rand + rand + rand + rand + rand + rand + rand - 6.0
		#return Math.rand_gauss # More accurate alternative
	end

	def generate_maxwell_velocity(atom_i, temperature)
		#: Generates random velocity for one atom from Maxwell distribution at given temperature
		atom = @geometry[atom_i]
		o = Coordinate.new
		3.times {|i|
			o[i] = atom[i] - rand_gauss * Math.sqrt(BOLTZMANN * temperature / atom.mass) * @dt
		}
		@oldcoord[atom_i] = o
		return nil
	end

	def generate_maxwell_velocities(init_temp)
		#: Generates random velocity for all atoms in the geometry from Maxwell distribution at given temperature
		@geometry.each_index {|i|
			generate_maxwell_velocity(i, init_temp)
		}
		# Save temperature
		t = temperature
		# Remove COM movement
		remove_com_translation
		remove_com_rotation
		# Get new temperature
		t_new = temperature
		# Rescale velocities to get the original temperature
		factor = Math.sqrt(t/t_new)
		rescale_velocities(factor)
	end

	#=======================================================================
	# Velocities manipulation
	#=======================================================================
	
	def rescale_velocities(factor)
		@geometry.each_with_index {|atom, i|
			@oldcoord[i].x = (1.0 - factor) * atom.x + factor * @oldcoord[i].x
			@oldcoord[i].y = (1.0 - factor) * atom.y + factor * @oldcoord[i].y
			@oldcoord[i].z = (1.0 - factor) * atom.z + factor * @oldcoord[i].z
		}
	end

	#=======================================================================
	# Thermodynamic variables
	#=======================================================================

	def kinetic_energy # => Float
		sum = 0.0
		@geometry.each_with_index {|atom, i|
			sum += atom.mass * ((atom.x-@oldcoord[i].x)**2 + (atom.y-@oldcoord[i].y)**2 + (atom.z-@oldcoord[i].z)**2)
		}
		return sum * 0.5 / @dt2
	end

	def temperature # => Float
		return kinetic_energy * 2.0 / (3.0 * @geometry.size * BOLTZMANN) 
	end

	def kinetic_and_temperature # => Array [kinetic_energy, temperature]
		ek = kinetic_energy
		t = ek * 2.0 / (3.0 * @geometry.size * BOLTZMANN)
		return [ek,t]
	end

	#=======================================================================
	# Center of mass motion removal
	#=======================================================================
	
	# The calculations here intentionaly don't use more advanced implementation
	# of algebra to save time.
	
	def com_shift # => [Coordinate, Coordinate]
		#: Returns centers of mass in this and previous step
		com     = Coordinate.new
		com_old = Coordinate.new
		totalmass = 0.0
		@geometry.each_with_index {|atom, i|
			totalmass += atom.mass
			com.plus!(atom * atom.mass)
			com_old.plus!(@oldcoord[i] * atom.mass)
		}
		com = com / totalmass
		com_old = com_old / totalmass
		return [com, com_old]
	end

	def remove_com_translation
		#: Translational velocity of the whole system is subtracted
		#: from atomic velocities.
		com, oldcom = com_shift
		shift = oldcom - com
		@geometry.each_with_index {|atom, i|
			@oldcoord[i].minus!(shift)
		}
		return nil
	end

	def angular_momentum(com) # => Coordinate
		#: Calculates angular momentum of the system, providing its center of mass
		#: is in origin of coordinates
		l = Coordinate.new
		@geometry.each_with_index{|atom, i|
			velocity = (atom - @oldcoord[i])/@dt
			l.plus!((atom - com).cross_product(velocity * atom.mass))
		}
		return l
	end

	def moment_of_inertia(com) # => [[ixx,ixy,ixz],[ixy,iyy,iyz],[ixz,iyz,izz]]
		#: Calculates moment of inertia of the system with respect to the
		#: supplied center of mass

		#: Returns array, not Matrix (for more efficiency)
		ixx = ixy = ixz = iyy = iyz = izz = 0.0
		@geometry.each { |atom|
			at = atom - com
			ixx += atom.mass * (at.y * at.y + at.z * at.z)
			iyy += atom.mass * (at.x * at.x + at.z * at.z)
			izz += atom.mass * (at.y * at.y + at.x * at.x)
			ixy -= atom.mass * at.x * at.y
			ixz -= atom.mass * at.x * at.z
			iyz -= atom.mass * at.z * at.y
		}
		return [[ixx,ixy,ixz],[ixy,iyy,iyz],[ixz,iyz,izz]]
	end

	def angular_velocity(com) # => Coordinate
		#: Calculates angular ivelocity of the system, providing its center of mass
		#: is in origin of coordinates
		jx,jy,jz = angular_momentum(com).to_a
		ixx,ixy,ixz,ixy,iyy,iyz,ixz,iyz,izz = moment_of_inertia(com).flatten
		# invert the moment of inertia -> "i1"
		det = ixx*iyy*izz + 2.0 * ixy*iyz*ixz - ixz*iyy*ixz - iyz*iyz*ixx - izz*ixy*ixy
		i1xx = (iyy*izz - iyz*iyz) / det
		i1yy = (ixx*izz - ixz*ixz) / det
		i1zz = (iyy*ixx - ixy*ixy) / det
		i1xy = (ixz*iyz - ixy*izz) / det
		i1xz = (ixy*iyz - ixz*iyy) / det
		i1yz = (ixz*ixy - ixx*iyz) / det
		# calculate center-of-mass angular velocity "w"
		wx = i1xx * jx + i1xy * jy + i1xz * jz
		wy = i1xy * jx + i1yy * jy + i1yz * jz
		wz = i1xz * jx + i1yz * jy + i1zz * jz
		return Coordinate.new(wx,wy,wz)
	end

	def remove_com_rotation # => nil
		#: Removes rotational motion of the system.
		com = @geometry.center_of_mass
		w = angular_velocity(com)
		@geometry.each_with_index { |atom, i|
			@oldcoord[i].plus!(w.cross_product(atom - com)*@dt)
		}
		return nil
	end

	#=======================================================================
	# Partial MD
	#=======================================================================
	
	def selection_gradient
		#: Returns gradient for the moving region
		if @atomlist
			return @gradient.for_atomlist(@atomlist)
		else
			return @gradient
		end
	end
	
	#=======================================================================
	# Integrators
	#=======================================================================

	def integrator_verlet
		#: Basic Verlet algorithm updating coordinates and velocities
		#: of all atoms in the system.
		selection_grad = selection_gradient
		@geometry.each_with_index {|atom, i|
			# x += x - xold - dt^2 * gradx / mass
			tempcoord = atom.to_coordinate
			atom.plus!(atom - @oldcoord[i] - selection_grad[i] * (@dt2 / atom.mass))
			@oldcoord[i] = tempcoord
		}
		return nil
	end

	#=======================================================================
	# History
	#=======================================================================

	def velocities_full
		if @atomlist
			velo_small = velocities
			velo_full = []
			zero = Coordinate.new(0,0,0)
			@geometry_full.each_index{|i| velo_full[i] = zero}
			@atomlist.each_index{|i|
				velo_full[atomlist[i]] = velo_small[i]
			}
			return velo_full
		else
			return velocities
		end
	end

	#=======================================================================
	# Variable temperature during simulation
	#=======================================================================
	def interpolate_temperature(cycle)
		differential = @settings[:temperature_target] - @settings[:temperature]
		point = (cycle + 1).to_f / @settings[:maxcycles]
		return @settings[:temperature] + differential * point
	end
end


