require "classes/internal_coordinates/z_matrix.rb"

module ActionsZmat
	def build_zmat
		# No geometry is written
		@write_geo = false

		# Build connectivity
		connectivity = Connectivity.new(@geometry)

		# Z-matrix 
		#zmat = ZMatrix.from_geometry_keep_order(@geometry)
		if @settings[:geometry_zmat_bondlist].size > 0
			bondlist = []
			@settings[:geometry_zmat_bondlist].each{|bond_sel|
				atomlist = @geometry.atomlist_from_selection(bond_sel)
				if atomlist.size != 2
					Cuby::error "Bond selection '#{bond_sel}' should select two atoms"
				end
				bondlist << atomlist
			}
			zmat = ZMatrix.from_geometry(@geometry, bondlist)
		else
			zmat = ZMatrix.from_geometry(@geometry)
		end

		case @settings[:geometry_zmat_format]
		when :gaussian
			puts zmat.to_s
		when :gaussian_vars
			zmat.build_variables!
			puts zmat.to_s
		when :mopac
			puts zmat.to_s_mopac(:user_bond)
		end
		puts
	end
end
