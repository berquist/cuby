require "protocols/scan/classes/scan_generator"

class ScanGeneratorGridXyz < ScanGenerator

	SETUP_KEYS = ["grid_x", "grid_y", "grid_z", "selection_move", "rvdw_min", "rvdw_max"]

	def initialize(settings)
		check_setup(settings, SETUP_KEYS)

		@ax = range_from_value(settings[:scan_generator_setup]["grid_x"], "grid_x")
		@ay = range_from_value(settings[:scan_generator_setup]["grid_y"], "grid_y")
		@az = range_from_value(settings[:scan_generator_setup]["grid_z"], "grid_z")
		@move_selection =  settings[:scan_generator_setup]["selection_move"].to_s
		@rvdw_min =  settings[:scan_generator_setup]["rvdw_min"].to_f
		@rvdw_max =  settings[:scan_generator_setup]["rvdw_max"].to_f
	end

	def iterate(geometry)
		# separate the geometry
		g_move = geometry.geometry_from_selection(@move_selection)
		l_move = geometry.atomlist_from_selection(@move_selection)
		l_static = geometry.atomlist_from_selection("%not(#{@move_selection})")

		# Centered copy of the moving geometry
		g_move_centered = g_move.deep_copy
		center = g_move.center
		g_move_centered.translate!(-center)

		count = 0
		@ax.each{|x| @ay.each{|y| @az.each{|z|
			crd = Coordinate[x,y,z]
			g_move.each_index{|i|
				g_move[i].set_coord(g_move_centered[i] + crd)
			}

			distance, ia, ib = geometry.closest_contact(l_move, l_static)
			rvdw = PeriodicTable.vdw_radius(geometry[ia].element) + PeriodicTable.vdw_radius(geometry[ib].element)
			if distance > rvdw * @rvdw_min && distance < rvdw * @rvdw_max
				str = sprintf("scan: %15.8f%15.8f%15.8f",x,y,z)
				yield geometry, str
				count += 1
			end
		}}}
		return count
	end

end
