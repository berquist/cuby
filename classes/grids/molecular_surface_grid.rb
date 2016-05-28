require "classes/grids/spherical_grid_spiral"

class MolecularSurfaceGrid
	attr_reader :coords

	def initialize(geometry, n, vdw_scale = 1.0)
		@coords = []

		# Precalculate scaled vdW radii
		radii = []
		geometry.each_index{|i|
			radii[i] = PeriodicTable.vdw_radius(geometry[i].element) * vdw_scale
		}

		geometry.each_index{|i|
			atom = geometry[i]
			atomgrid = SphericalGridSpiral.new(radii[i], n)
			atomgrid.iterate_xyz(atom){|crd|
				save = true
				geometry.each_index{|j|
					next if i == j
					if crd.distance(geometry[j]) < radii[j]
						save = false
						break
					end
				}
				@coords << crd if save
			}
		}
	end
end
