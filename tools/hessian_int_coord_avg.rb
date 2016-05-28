#!/usr/bin/ruby
$:.unshift(File.dirname(__FILE__)+"/..").uniq! # Magic line to add path to the executable to $RUBYLIB
require "classes/cuby.rb" # Load Cuby classes and libraries

require "yaml"
require "classes/objects/hessian_internal_coordinates.rb"

if ARGV[0]
	read_hic = HessianInternalCoordinates.from_yaml(ARGV[0])
	write_hic = HessianInternalCoordinates.new
	
	read_hic.diag_elem_array.each do |read_item|
		incl = write_hic.include_by_type(read_item)
		if incl[0]
			write_item = write_hic.diag_elem_array[incl[1]]
			write_hic.diag_elem_array[incl[1]].hess_value = (write_item.hess_value * write_item.weight + read_item.hess_value * read_item.weight) / (write_item.weight + read_item.weight)
			write_hic.diag_elem_array[incl[1]].coord_value = (write_item.coord_value * write_item.weight + read_item.coord_value * read_item.weight) / (write_item.weight + read_item.weight)
			write_hic.diag_elem_array[incl[1]].weight = write_item.weight + read_item.weight
		else
			write_hic.diag_elem_array << read_item
		end
	end
	
	puts write_hic.diag_elem_array.to_yaml
	
else
	Cuby::error("You must provide an argument (input filename).")
end
