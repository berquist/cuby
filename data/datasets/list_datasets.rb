#!/usr/bin/env ruby
$:.unshift(File.dirname(`which cuby4`)).uniq! # Magic line to add path to the executable to $RUBYLIB
require "yaml"
require "protocols/dataset/classes/data_set"

dataset_path = File.dirname(`which cuby4`)+"/data/datasets"

count = 1
Dir.entries(dataset_path).sort.each{|fn|
	next unless fn =~ /\.yaml$/
	f = File.open(dataset_path+"/"+fn)
	set = YAML.load(f)
	f.close

	name = set.description.name

	line = "<tr>"
	line += "<td><strong>" + name + "</strong></td>"
	line += "<td>" + set.description.text
	set.description.references.each_pair{|doi, ref|
		line += "<a href=\"http://dx.doi.org/#{doi}\" title=\"#{ref}\">[#{count}]</a>"
		count += 1
	}
	line += "</td>"
	puts line
}


