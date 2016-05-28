
module ProtocolShellScript
	#=======================================================================
	# Protocol header
	#=======================================================================
	PROTOCOL_TYPE = :tools
	PROTOCOL_DESCRIPTION = "Shell commands execution"
	#=======================================================================

	def print_header_protocol
		# Print header
		Cuby.log.print_cuby_logo(:comment => false, :text => 'Shell script')
	end

	def prepare_protocol
	end

	def run_protocol(queue)
		system(@settings[:shell_commands])
	end

	def print_protocol
	end

	def cleanup_protocol
	end

	def results
		return nil
	end

end
