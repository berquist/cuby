#!/usr/bin/env ruby
$:.unshift(File.absolute_path(Dir.pwd + "/../../..")).uniq! # Magic line to add path to the executable to $RUBYLIB
require "classes/cuby.rb" # Load Cuby classes and libraries
require "classes/core/cuby_loader.rb"
require "classes/grids/molecular_surface_grid.rb"

g = Geometry.from_file(ARGV[0])

grid = MolecularSurfaceGrid.new(g, 200)

grid.coords.each{|crd|
	g << Atom.from_coordinate(:X, crd)
}

g.write_xyz
