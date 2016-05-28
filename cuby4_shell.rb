#!/usr/bin/env ruby
$:.unshift(File.dirname(__FILE__)).uniq! # Magic line to add path to teh executable to $RUBYLIB
require "classes/cuby.rb" # Load Cuby classes and libraries

# Debug messages
if ARGV[0] == '-d'
	Cuby.log.logs[0].verbosity = :debug 
	ARGV.shift
end

require "irb"

# Load the farmework
Cuby.load_framework
# Configure cuby using default settings
Cuby.configure(Settings.new)

Cuby.log.print_cuby_logo(:text => "Interactive shell")

if __FILE__ == $0
  IRB.start(__FILE__)
else
  # check -e option
  if /^-e$/ =~ $0
    IRB.start(__FILE__)
  else
    IRB.setup(__FILE__)
  end
end
