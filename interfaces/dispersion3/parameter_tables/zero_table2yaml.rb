# This scripts generates the yaml data for ../dispersion3_parameters.yaml
# from the table of parameters available at Stefan Grimme's website
# http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=DFT-D3&lang=english

filename = "grimme_zero_damping.txt"
indent = "              "

IO.readlines(filename).each{|line|
	func, s6, sr6, s8 = line.strip.split
	func = func.downcase.gsub("-","")

	puts indent + func + ":"
	puts indent + "  d3_s6: " + s6 unless s6.to_f == 1.0
	puts indent + "  d3_sr6: " + sr6
	puts indent + "  d3_s8: " + s8

}
