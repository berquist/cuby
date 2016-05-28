#!/usr/bin/env ruby

# Print some help if no argument provided
if ARGV.size == 0
puts "
Usage:
dataset_from_table.rb table.dat > set.yaml

The table should contain following fields, separated by tabulators:
number, name, shortname, symmetry, energy, group, tags, geometry
\n"
exit
end

TEMPLATE = "- !ruby/object:ProtocolDataset::DataSetItem
  name: %name%
  shortname: %shortname%
  geometry: %geometry%
  reference_value: %energy%
  setup: {%setup%}
  group: %group%
  tags: \"%tags%\""


table = IO.readlines(ARGV[0])

table.each{|line|
	number, name, shortname, symmetry, energy, group, tags, geometry, sel_a, sel_b, c_a, c_b = line.strip.split("\t")
	next if number == "number"

	number = "%02d" % number.to_i

	s = TEMPLATE.dup
	s.gsub!("%name%", "#{number} #{name}")
	s.gsub!("%shortname%", "#{number}_#{shortname}")
	if geometry && geometry.strip != ""
		if ARGV[1]
			s.gsub!("%geometry%", "#{ARGV[1]}:#{geometry}")
		else
			s.gsub!("%geometry%", "#{geometry}")
		end
	else
		s.gsub!("%geometry%", "#{ARGV[1]}:#{number}")
	end
	
	c_a = 0 unless c_a
	c_b = 0 unless c_b

	if sel_a
		if sel_a =~ /^%/ || sel_b =~ /^%/
			sel_a = '"'+sel_a+'"'
			sel_b = '"'+sel_b+'"'
		end

		setup = "
    molecule_a:
      selection: #{sel_a}
      charge: #{c_a}
    molecule_b:
      selection: #{sel_b}
      charge: #{c_b}"
		s.gsub!("{%setup%}", "#{setup}")
	end
	s.gsub!("%energy%", energy)
	s.gsub!("%group%", group)
	s.gsub!("%tags%", tags)

	puts s
	puts
}
