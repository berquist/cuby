#!/usr/bin/env ruby
$:.unshift(File.dirname(__FILE__)).uniq! # Magic line to add path to the executable to $RUBYLIB
require "classes/cuby.rb" # Load Cuby classes and libraries
require 'drb'
require 'socket'

require "classes/calculation/calculation_external.rb"

# Load the farmework
Cuby.load_framework
# Configure cuby using default settings
Cuby.configure(Settings.new)

class CubyServer
	#attr_reader :calculation

	def initialize
		@workdir = Dir.pwd
		@hostname = Socket.gethostname
	end

	def workdir
		return @workdir
	end

	def hostname
		return @hostname
	end

	def calculate(calculation, interface_file, working_directory)
		require(interface_file)
		calculation = Marshal.load(calculation)
		Dir.chdir(working_directory)
		puts "Running calculation #{calculation.name}"
		puts "in #{working_directory}"
		results = calculation.calculate_interface
		puts "Finished"
		puts
		$stdout.flush
		return results
	end

	def stop_server
		puts "Remote request to exit"
		puts
		$stdout.flush
		# New thread is created to delay the shutdown so that the
		# communication can be finished
		Thread.new{ 
			sleep(5)
			exit
		}
	end
end

if ARGV.size > 0
	port = ARGV[0]
	DRb.start_service("druby://:#{port}", CubyServer.new)
else
	DRb.start_service(nil, CubyServer.new)
end


# Print URI, this is passed to the client
puts "Server URI:"
puts DRb.uri
puts
$stdout.flush

# Print the uri into a separate file
File.open("server_uri","w+"){|f|
	        f.puts DRb.uri
}


DRb.thread.join
