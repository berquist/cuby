################################################################################
# 
# CUBY libraries loader & module Cuby
#
################################################################################

# This file provides the wapper for cuby applications and loads the objects used
# by cuby: Log files for output, error handling, timer, global configuration etc.

require "classes/core/logfiles.rb"
require "fileutils.rb"

class CubyError < RuntimeError
	#: This exception class is used when the program encounters wrong input from
	#: user or some other well-defined failure occurs. 
end

module Cuby
	#=======================================================================
	# Environment setup
	#=======================================================================
	
	def Cuby.load_framework(options = {})
		#: Load all cuby classes and libraries and initialize Settings class.
		#: Called from the code so that Cuby module can be configured
		#: before classes are loaded.
		require "classes/core/cuby_loader.rb"

		# Load keyword list from file
		Settings.load_keyword_list(Cuby.install_dir + "/classes/keywords.yaml")
		Settings.load_keyword_list_interfaces(Cuby.install_dir + "/interfaces")
		Settings.load_keyword_list_protocols(Cuby.install_dir + "/protocols")

		# Cuby config files works by redefining default values of Settings
		# keywords
	
		# 1) Global config - located one level above the cuby installation
		#    (to prevent overwriting it on update)
		configfiles = 0
		unless options[:skip_global_config]
			config_fn = File.expand_path(Cuby.install_dir + "/../cuby4_config.yaml")
			if  FileTest.exist?(config_fn)
				configfiles += 1 
				Settings.load_keyword_defaults(config_fn)
			end
		end

		# 2) User config file
		unless options[:skip_user_config]
			local_config_fn = File.expand_path('~/.cuby4')
			if FileTest.exist?(local_config_fn)
				configfiles +=  1
				Settings.load_keyword_defaults(local_config_fn)
			end
		end

		# Warn when there are no config files at all
		if configfiles == 0 && !options[:skip_global_config] && !options[:skip_user_config]
			$stderr.puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			$stderr.puts ""
			$stderr.puts " WARNING - no configuration file found"
			$stderr.puts ""
			$stderr.puts " Global config file should be located at #{File.expand_path(Cuby.install_dir + "/../cuby4_config.yaml")}"
			$stderr.puts " User-specific config is read from ~/.cuby4"
			$stderr.puts ""
			$stderr.puts " Cuby will run but all the interfaces have to be configured from the input file."
			$stderr.puts ""
			$stderr.puts "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		end
	end

	#=======================================================================
	# Cuby application wrapper
	#=======================================================================

	def Cuby.application
		#: Wrapper for all Cuby executables. Provides nice exception handling
		#: and run time counter to the block it runs.
		begin
			Cuby.start_timer
			yield
			Cuby.stop_timer
			Cuby.print_timer if @@print_timer
		rescue
			if $!.class == CubyError
				Cuby.error_handle
			else
				Cuby.exception_handle
			end
		end
	end

	#=======================================================================
	# Log file	
	#=======================================================================

	# Create log file with default settings
	@@log = Logfiles.new
	@@log.add_log(:normal, :stdout, :stderr)

	def Cuby.log
		#: Accessor to the log object
		return @@log
	end

	#=======================================================================
	# Calculation queue
	#=======================================================================
	
	@@queue = nil

	def Cuby.queue
		#: Accessor for the queue
		return @@queue
	end
	
	#=======================================================================
	# Config
	#=======================================================================
	
	def Cuby.configure(settings)
		#: Reads global settings stored in Cuby module from
		#: a settings object.

		# Set verbosity
		@@log.logs[0].verbosity = settings[:verbosity]

		# Timer printing
		@@print_timer = settings[:print].include?(:timing) && ! settings[:queue_submit]

		# Create default queue
		@@queue = CalculationQueue.new(settings[:cuby_threads])

		# Initialize remote calculation execution if needed
		CalculationRemote::init_remote(settings)

		# Recommendation exception handling
		@@ignore_recommendations = settings[:ignore_recommendations]

		# Local data directory
		Cuby.local_data_init(settings[:cuby_local_data_dir])
	end

	#=======================================================================
	# Finalization
	#=======================================================================
	
	def Cuby.finalize
		CalculationRemote.stop_servers
	end

	#=======================================================================
	# Installation directory
	#=======================================================================

	# Installation directory is determined from the path of the executable
	@@install_dir = File.dirname(__FILE__).gsub(/\/classes$/,"")

	def Cuby.install_dir
		#: Getter for install_dir
		return @@install_dir
	end

	#=======================================================================
	# Local data directory
	#=======================================================================
	
	def Cuby.local_data_init(dirname)
		dirname = File.expand_path(dirname)

		if FileTest.directory?(dirname)
			# Directory exists
			@@local_data_dir = dirname
		elsif FileTest.exist?(dirname)
			# File exists, but is not a directory
			Cuby::error "Keyword 'cuby_local_data_dir' points to a file, not a directory as it is supposed to (#{dirname})."
		else
			# Create new one
			 @@local_data_dir = dirname
			Cuby::log.puts_debug "Creating local data directory #{@@local_data_dir}"
			begin
				FileUtils.mkdir(dirname)
			rescue
				Cuby::error "Failed to create directory specified by keyword 'cuby_local_data_dir' (#{dirname})."
			end
		end
	end

	def Cuby.local_data_dir
		#: Getter for local_data_dir
		return @@local_data_dir
	end

	#=======================================================================
	# Error handling
	#=======================================================================
	
	# Initialize storage of warnings
	@@warnings_printed = []
	
	def Cuby.error(message, context = nil)
		#: Synonym to raise CubyError
		@@cuby_error_context = context
		raise(CubyError, message)
	end

	def Cuby.warning(error_string, repeat = :always)
		#: Print warning
		if repeat == :once
			hash = error_string.crypt('cuby')
			return false if @@warnings_printed.include?(hash)
			@@warnings_printed << hash
		end

		$stderr.print "WARNING: "
		line = 0
		error_string.each_line{|s|
			if line == 0
				$stderr.puts "#{s}"
			else
				$stderr.puts "         #{s}"
			end
			line += 1
		}
		$stderr.flush
		return true
	end

	def Cuby.recommendation(message)
		#: Method for nice printing of CubyRecommendationError exceptions used by Cuby.application
		@@log.flush
		@@log.pute "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		if @@ignore_recommendations
			@@log.pute "! RECOMMENDATION (IGNORED)"
		else
			@@log.pute "! RECOMMENDATION: "
		end
		message.each_line{|s| @@log.pute "! #{s}" }
		if @@ignore_recommendations
			@@log.pute "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
			@@log.pute ""
			@@log.flush
			return
		end
		@@log.pute "!"
		@@log.pute "! Cuby thinks you are attempting to do something technically possible but wrong."
		@@log.pute "! If you are sure what you are doing, restart the calculation with keyword"
		@@log.pute "! ignore_recommendations: yes"
		@@log.pute "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@@log.pute ""
		@@log.flush
		exit(1)
	end


	def Cuby.exception_handle(details = $!)
		#: Method for nice printing of exceptions used by Cuby.application
		@@log.flush
		@@log.pute "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@@log.pute "! ERROR (#{details.class}): "
		unless details.message.nil?
			details.message.each_line{|s| @@log.pute "! #{s}"}
		end
		@@log.pute "!"
		@@log.pute "! Backtrace:"
		details.backtrace.each{|s| @@log.pute "! #{s}"}
		@@log.pute "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@@log.pute ""
		@@log.flush
		exit(1)
	end

	def Cuby.error_handle(details = $!)
		#: Method for nice printing of CubyError exceptions used by Cuby.application
		@@log.flush
		@@log.pute "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		case @@cuby_error_context
		when Calculation
			@@log.pute "! ERROR in calculation #{@@cuby_error_context.name}, interface #{@@cuby_error_context.settings[:interface]}: "
		else
			@@log.pute "! ERROR: "
		end
		unless details.message.nil?
			details.message.each_line{|s| @@log.pute "! #{s}" }
		end
		@@log.pute "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
		@@log.pute ""
		@@log.flush
		exit(1)
	end

	#=======================================================================
	# Timer
	#=======================================================================

	# The execution time of the Cuby::application block is measured by
	# the following timer
	
	@@print_timer = false

	def Cuby.start_timer
		#: Start Cuby module timer
		@@start_time = Time.now
	end

	def Cuby.stop_timer
		#: Stop Cuby module timer
		@@stop_time = Time.now
	end

	def Cuby.print_timer
		# Print timing to the log file
		@@log.puts
		@@log.puts "Timing (in various units)"
		@@log.puts(sprintf("   Seconds: %15.2f\n",(@@stop_time - @@start_time).to_f))
		@@log.puts(sprintf("   Minutes: %15.2f\n",(@@stop_time - @@start_time).to_f/60))
		@@log.puts(sprintf("   Hours:   %15.2f\n",(@@stop_time - @@start_time).to_f/60/60))
		@@log.puts(sprintf("   Days:    %15.2f\n",(@@stop_time - @@start_time).to_f/60/60/24))
	end

end
