#!/usr/bin/env ruby
$:.unshift('/home/rezac/cuby4')
require "classes/cuby.rb" # Load Cuby classes and libraries
require "classes/core/cuby_loader.rb"
require "classes/geometry/hybridization.rb"

g = Geometry.from_file(ARGV[0])

hyb = Hybridization.new(g)
hyb.evaluate_grimme(:PA2009, true)

puts "hybridization"
puts hyb
puts

puts	"numerical derivatives"
d = 0.0001
derivatives = []
g.each_index{|i|
	derivatives[i] = Coordinate.new
	3.times{|c|
		g[i][c] += d
		hp =  hyb.grimme_atom(i,:PA2009, true)
		g[i][c] -= d * 2
		hm = hyb.grimme_atom(i,:PA2009, true)
		g[i][c] += d # restore
		derivatives[i][c] = (hp - hm) / d / 2
	}
}
derivatives.each_index{|i|
	d = derivatives[i]
	puts "#{i}:\t#{'%.8f'%d.x}\t#{'%.8f'%d.y}\t#{'%.8f'%d.z}"
}

g = Geometry.from_file(ARGV[0])

hyb = Hybridization.new(g)
puts
puts "analytical derivatives"
hyb.evaluate_grimme(:PA2009, true, true)
derivatives_a = hyb.gradient
derivatives_a.each_index{|i|
	d = derivatives_a[i]
	puts "#{i}:\t#{'%.8f'%d.x}\t#{'%.8f'%d.y}\t#{'%.8f'%d.z}"
}

puts
puts "difference"
derivatives_a.each_index{|i|
	d = derivatives_a[i] - derivatives[i]
	puts "#{i}:\t#{'%.8f'%d.x}\t#{'%.8f'%d.y}\t#{'%.8f'%d.z}"
}
