class GeoGridCubic
	# The geometry is split into cubes of defined size, each cube is an aray of
	# indexes pointing to the original geometry

	attr_reader :size_x, :size_y, :size_z
	attr_reader :cubes
	
	def initialize(geometry, cube_size = 5)
		# Get the max/min coordinates
		coord_max = Coordinate[-Float::MAX, -Float::MAX, -Float::MAX]
		coord_min = Coordinate[Float::MAX, Float::MAX, Float::MAX]
		geometry.each{|atom|
			coord_max.x = atom.x if atom.x > coord_max.x
			coord_max.y = atom.y if atom.y > coord_max.y
			coord_max.z = atom.z if atom.z > coord_max.z

			coord_min.x =  atom.x if atom.x < coord_min.x
			coord_min.y =  atom.y if atom.y < coord_min.y
			coord_min.z =  atom.z if atom.z < coord_min.z
		}

		# Get grid size
		boxsize = coord_max - coord_min
		@size_x = (boxsize.x / cube_size).ceil
		@size_y = (boxsize.y / cube_size).ceil
		@size_z = (boxsize.z / cube_size).ceil

		# Initialize the cubes
		@cubes = []
		@size_x.times{|ix|
			@cubes[ix] = []
			@size_y.times{|iy|
				@cubes[ix][iy] = []
				@size_z.times{|iz|
					@cubes[ix][iy][iz] = []
				}
			}
		}

		# Fill the cubes
		geometry.each_index{|i|
			# Get cube index from atom coordinate
			dc = geometry.at(i) - coord_min
			ix = (dc.x / cube_size).floor
			iy = (dc.y / cube_size).floor
			iz = (dc.z / cube_size).floor
			# Save atom
			@cubes[ix][iy][iz] << i
		}
	end

	def each_index
		@size_x.times{|ix|
			@size_y.times{|iy|
				@size_z.times{|iz|
					yield ix, iy, iz
				}
			}
		}
	end

	def each
		each_index{|ix, iy, iz| yield @cubes[ix][iy][iz]}
	end

	ADJACENT_DIRECTIONS = [[1,0,0], [1,1,0], [0,1,0], [-1,1,0], 
		[0,0,1], [1,0,1], [1,1,1], [0,1,1], [-1,0,1], [-1,-1,1], [0,-1,1], [-1,1,1], [1,-1,1]]

	def each_adjacent_pair_indexes
		each_index{|ix, iy, iz|
			ADJACENT_DIRECTIONS.each{|direction|
				dx, dy, dz = direction
				nx = ix + dx
				next if nx < 0
				next if nx >= @size_x
				ny = iy + dy
				next if ny < 0
				next if ny >= @size_y
				nz = iz + dz
				next if nz < 0
				next if nz >= @size_z
				yield ix, iy, iz, nx, ny, nz
			}
		}
	end

	def each_adjacent_pair(exclude_self = false)
		each_adjacent_pair_indexes{|ix, iy, iz, nx, ny, nz|
			yield @cubes[ix][iy][iz], @cubes[nx][ny][nz]
		}
	end
end

