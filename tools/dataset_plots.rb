
$:.unshift(File.dirname(__FILE__)+"/..").uniq! # Magic line to add path to the executable to $RUBYLIB
require "yaml"
require "/home/rezac/cuby4/protocols/dataset/classes/data_set.rb"

# Help
if ARGV.size != 3
	puts "Use:"
	puts ""
	puts "dataset_plots.rb dataset.yaml scan_start scan_end"
	puts ""
	puts "where scan_start and scan_end are numbers identifying the first and last point of the scan,"
	puts "these numbers should be the last word of the item name"
	puts ""
	exit 1
end

# print header
puts "#==============================================================================="
puts "# Plots"
puts "#==============================================================================="
puts ""
puts "plots:"

# load dataset file
set = YAML.load(File.open(ARGV[0]))

first = ARGV[1].to_f
last = ARGV[2].to_f

i_f = nil
xes = nil
set.items.each_index{|i|
	item = set.items[i]
	name = item.name.strip.split(/ +|_+/)
	if name.last.to_f == first
		i_f = i
		xes = []
	end
	xes << name.last.to_f
	if name.last.to_f == last
		longname = item.name.split
		longname.pop
		longname = longname.join(" ")
		shortname = item.shortname.split("_")
		shortname.pop
		shortname = shortname.join("_")
		puts "- !ruby/object:ProtocolDataset::DataSetPlot"
		puts "  name: \"#{longname}\""
		puts "  filename: #{shortname}"
		puts "  first: #{i_f}"
		puts "  last: #{i}"
		puts "  x_values: [#{xes.join(', ')}]"
	end
}

#--------------------------------------------------
# - !ruby/object:ProtocolDataset::DataSetPlot
#   name: "01 Water ... Water"
#   filename: 01_Water-Water
#   first: 0
#   last: 7
#   x_values: [0.9, 0.95, 1.0, 1.05, 1.1, 1.25, 1.5, 2.0]
# 
#-------------------------------------------------- 
