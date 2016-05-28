
class Hessian < Algebra::Matrix

	#=======================================================================
	# Mass weighting matrix
	#=======================================================================

	def Hessian.mass_weight_matrix(geometry)
		m = Algebra::Matrix.zero(geometry.size*3)
		(geometry.size*3).times {|i|
		        m[i,i] = 1.0 / (geometry[i/3].mass ** 0.5) # modify the diagonal element
		}
		return m
	end

	def Hessian.mass_unweight_matrix(geometry)
		m = Algebra::Matrix.zero(geometry.size*3)
		(geometry.size*3).times {|i|
		        m[i,i] = (geometry[i/3].mass ** 0.5) # modify the diagonal element
		}
		return m
	end

	#=======================================================================
	# Projector matrices
	#=======================================================================

	def Hessian.transrot_removal_projector(geometry, mass_weighted = false, remove_translations = true)
		#: Builds projector matrix for removal of translation and rotation.
		
		Cuby::error "Translation/rotation projector can not be calculated for empty geometry" if geometry.size == 0

		# move center to [0,0,0]
		if mass_weighted
			com = geometry.weighted_center(:mass)
		else
			com = geometry.weighted_center(:unit)
		end
		geometry.translate!(-com)

		# Array of Vectors to remove
		remove_vectors = []

		# Linear or other?
		axes = []
		eigvectors, eigenvalues = geometry.moment_of_inertia_diagonalized
		if eigenvalues[0].abs < 1.0e-7
			# Linear molecule
			# The two roation exes perpendicular to the molecular axis [eigvectors[0]:
			axes << Coordinate[eigvectors[1][0], eigvectors[1][1], eigvectors[1][2]]
			axes << Coordinate[eigvectors[2][0], eigvectors[2][1], eigvectors[2][2]]
		else
			# General case: use cartesian coordinates
			3.times{|a|
				axis = Coordinate[0,0,0]
				axis[a] = 1.0
				axes << axis
			}
		end

		# Rotations
		axes.each{|axis|
			vector = Vector.zero(geometry.size * 3)
			geometry.size.times {|i|
				c = geometry[i]
				if mass_weighted
					c *= geometry[i].mass ** 0.5
				end

				cross = c.cross_product(axis)
				vector[3*i]   = cross.x
				vector[3*i+1] = cross.y
				vector[3*i+2] = cross.z
			}
			remove_vectors << vector
		}

		if remove_translations
			# Add translations (prepend, so that thay are indexed as 0,1,2)
			3.times {|i| remove_vectors.unshift(Vector.zero(geometry.size * 3))}

			geometry.size.times {|i|
				# Translations
				3.times {|c|
					if mass_weighted
						remove_vectors[c][3*i + c] = geometry[i].mass ** 0.5
					else
						remove_vectors[c][3*i + c] = 1.0
					end
				}
			}
		end

		# Orthonormalize vectors
		Vector.orthonormalization_gram_schmidt(remove_vectors)

		# Construct projector
		c = Algebra::Matrix.columns(remove_vectors.map{|v| v.to_a})
		p = Algebra::Matrix.identity(geometry.size * 3) - c * c.transpose

		geometry.translate!(com)
		return p
	end

	def Hessian.rot_removal_projector(geometry, mass_weighted = false)
		return  Hessian.transrot_removal_projector(geometry, mass_weighted, false)
	end

	#=======================================================================
	# Transformation from calcultion geometry
	#=======================================================================

	def Hessian.from_calculation_coordinates(geo_cuby, geo_new, hess_new)
		# Finds an orthogonal transformation matrix (not necessarily only a rotation) that maps
		# the calculation geometry onto the input one, and applies this transformation to
		# the Hessian. The solution uses:
		# http://en.wikipedia.org/wiki/Orthogonal_Procrustes_problem

		com_i = geo_cuby.center_of_mass
		com_o = geo_new.center_of_mass

		gi = geo_cuby.translate(-com_i)
		go = geo_new.translate(-com_o)

		p = gi.to_matrix
		q = go.to_matrix
		a = p.transpose * q
		u, s, v = a.svd

		# 3x3 transformation matrix
		rotm = v.transpose * u.transpose

		# rotator matrix in 3n
		#--------------------------------------------------
		# rot2 = Hessian.zero(hess_new.m, hess_new.n)
		# geo_cuby.size.times{|i|
		# 	rot2.paste!(i*3,i*3,rotm)
		# }
		# return rot2.transpose * hess_new * rot2
		#-------------------------------------------------- 

		# Avoid construction of large matrix:
		# multiply each 3x3 block o hessian
		hess_rotated = Hessian.zero(hess_new.m, hess_new.n)
		rotm_t = rotm.transpose
		geo_cuby.size.times{|i|
			geo_cuby.size.times{|j|
				hess_rotated.paste!(i*3, j*3, rotm_t * hess_new.submatrix(i*3, j*3, 3, 3) * rotm)
			}
		}
		return hess_rotated
	end

	#=======================================================================
	# YAML support
	#=======================================================================
	
	# Copied from Matrix, class name changed
	
	require "yaml"

	begin
		yaml_engine = "psych" if YAML == Psych
		yaml_engine = "syck" if YAML == Syck
	rescue
		# ruby1.8 uses syck by default
		yaml_engine = 'syck'
	end

	if yaml_engine == 'syck'
		# Syck yaml engine

		yaml_as "tag:cuby.molecular.cz,2009:#{self}"

		def to_yaml(opts = {})
			YAML::quick_emit( self.object_id, opts ) do |out|
				out.map(taguri) do |map|
					arr = self.to_a
					arr.each{|a| def a.to_yaml_style; :inline; end}
					map.add('elements', arr)
				end
			end
		end

		def self.yaml_new(klass, tag, val)
			array = val['elements']
			return Hessian.from_array(array)
		end

	elsif yaml_engine == 'psych'
		# Psych yaml engine
	
		def encode_with(coder)
			coder.tag = "!cuby.molecular.cz,2009/Hessian"
			arr = self.to_a
			coder.style = Psych::Nodes::Mapping::BLOCK # BLOCK or FLOW
			#!# Would be nice to have rows styled as FLOW, but this will be more complicated
			coder['elements'] = arr
		end

		def init_with(coder)
			puts "init_with"
		end

		YAML::add_domain_type('cuby.molecular.cz,2009', 'Hessian') { |type, value|
			Hessian.from_array(value['elements'])
		}

	else
		raise "YAML engine can not be determined"
	end
end
