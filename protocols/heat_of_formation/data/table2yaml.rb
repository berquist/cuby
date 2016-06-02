#!/usr/bin/env ruby

data = IO.readlines(ARGV[0])
data.map!{|l| l.strip.split}

data.each{|line|
	ele, h0, h298, s298 = line
	h0 = "%.3f" % h0
	h298 = "%.3f" % h298 if h298 =~ /[0-9]/
	s298 = "%.3f" % s298

	puts ":#{ele}: !ruby/object:AtomThermochemistry
  h0: #{h0}
  h298_corr: #{h298}
  s0_298: #{s298}

"
}
