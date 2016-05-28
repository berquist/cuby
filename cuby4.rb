#!/usr/bin/env ruby
$:.unshift(File.dirname(__FILE__)).uniq! # Magic line to add path to the executable to $RUBYLIB
require "classes/cuby.rb" # Load Cuby classes and libraries

Cuby.application {
	# Load Cuby classes and libraries
	Cuby.load_framework

	# Parse commandline and build settings
	require_input = ARGV.size == 0
	settings = Settings.from_commandline(ARGV, require_input)

	# Configure Cuby module using the settings
	Cuby.configure(settings)

	if settings[:queue_submit]
		# Submit job to queue
		require "classes/queue_system/queue_system.rb"
		queue_system = QueueSystem.new(settings, ARGV)
		queue_system.submit!
	else
		# Run job
		job = Job.new("job", settings)
		job.print_header
		job.prepare
		job.run
		job.print
		job.cleanup

		# Write parsed input if requested
		settings.save_to_file("parsed_input.yaml") if settings[:write_parsed_input]

		# Print empty line as a separator in more verbose modes
		Cuby.log.puts_v(:normal,"")
	end

	Cuby.finalize # Cleanup etc
}
