module CalculationRemote

	#=======================================================================
	# Setup for remote calculations
	#=======================================================================
	
	def CalculationRemote.init_remote(settings)
		@@remote = settings[:remote_calculations]
		case @@remote
		when :no
			# nothing to be done
		when :drb
			CalculationRemote.init_drb(settings)
		end

		# Workdir handling
		@@remote_workdir_mode = settings[:remote_file_transfer] unless @@remote == :no

		# Stop servers when finished
		@@remote_stop_servers = settings[:remote_servers_close]
	end

	def CalculationRemote.init_drb(settings)
		# Load DRB
		require 'drb'

		# Initialize DRB
		DRb.start_service

		# Number of servers must be equal to number of cuby_threads
		unless settings[:remote_drb_servers].size == settings[:cuby_threads]
			Cuby::error("Number of remote servers must be equal to the number of parallel threads")
		end

		# Connect to the servers
		@@remote_drb_servers = []
		settings[:remote_drb_servers].each{|uri|
			begin
				calc_server = DRbObject.new(nil, uri)
				@@remote_drb_servers << calc_server
			rescue
				Cuby::error("Failed to connect to DRB server #{uri}")
			end
		}
	end

	def CalculationRemote.stop_servers
		return unless @@remote_stop_servers
		case @@remote
		when :drb
			@@remote_drb_servers.each{|server|
				server.stop_server
			}
		end
	end

	#=======================================================================
	# Working directory handling
	#=======================================================================
	
	def remote_copy_workdir_to_server(hostname, server_workdir)
		case @@remote_workdir_mode
		when :shared_wd
			# Return path to the current directory, should be identical on the server
			return Dir.pwd
		when :scp
			# Copy calculation dir to the server
			`scp -r #{calc_dir} #{hostname}:#{server_workdir}`
			# Return path to the server's working directory
			return server_workdir
		end
	end

	def remote_copy_workdir_from_server(hostname, server_workdir)
		case @@remote_workdir_mode
		when :shared_wd
			# Nothing to be done
		when :scp
			# Copy calculation dir from the server
			`scp -r #{hostname}:#{server_workdir}/#{calc_dir} .`
			# Remove the remote copy
			`echo 'rm -rf #{server_workdir}/#{calc_dir}' | ssh #{hostname} 2> /dev/null`
		end
	end
end
