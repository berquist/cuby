require "classes/core/logfile.rb"

class Logfiles
	attr_reader :logs

	def initialize
		@logs = []
	end

	def add_log(verbosity, file, errorfile)
		@logs << Logfile.new(verbosity, file, errorfile)
	end


	def puts(object = '')
		@logs.each{|log|
			log.puts(object)
		}
	end	

	def pute(object = '')
		@logs.each{|log|
			log.pute(object)
		}
	end	

	def flush
		@logs.each{|log|
			log.flush
		}
	end

	def puts_v(verbosity, object = '')
		@logs.each{|log|
			log.puts_v(verbosity, object)
		}
	end

	def puts_v_only(verbosity, object = '')
		@logs.each{|log|
			log.puts_v_only(verbosity, object)
		}
	end

	def puts_debug(object = '')
		@logs.each{|log|
			log.puts_debug(object)
		}
	end
	
	def print_cuby_logo(options = {})
		@logs.each{|log|
			log.print_cuby_logo(options)
		}
	end

end
