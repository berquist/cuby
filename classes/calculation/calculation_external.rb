# This module provides methods shared by interfaces calling external programs
# for the calculation

require 'fileutils'
require 'digest'

module CalculationExternal

	#=======================================================================
	# Interface checks
	#=======================================================================

	def check_method_availability
		unless @interface_module::METHODS.has_key?(@settings[:method])
				Cuby::error("Calculation \"#{@name}\":\nInterface #{@settings[:interface]} does not support method #{@settings[:method]}")
		end
	end

	def check_method_capabilities
		@what.each{|what|
			unless @interface_module::METHODS[@settings[:method]].include?(what)
				Cuby::error("Calculation \"#{@name}\":\nInterface #{@settings[:interface]} does not support calculation of #{what} by method #{@settings[:method]}")
			end
		}
		if @point_charges && !@interface_module::METHODS[@settings[:method]].include?(:point_charges)
			Cuby::error("Calculation \"#{@name}\":\nInterface #{@settings[:interface]} does not support calculation with extrenal point charges by method #{@settings[:method]}")
		end
	end

	#=======================================================================
	# Calculation directory
	#=======================================================================
	
	def calc_dir_short?
		# Switch between short or long notation
		return false
	end
	
	def calc_dir
		if calc_dir_short?
			return "j_" + Digest::MD5.hexdigest(@name)[0..13]
		else
			return @name + "_" + @interface_module::DIRECTORY
		end
	end

	def in_calc_dir(filename)
		return calc_dir + "/" + filename
	end

	def calc_dir_mkdir(input_file = nil, output_file = nil)
		if FileTest.exist?(calc_dir)
			if FileTest.directory?(calc_dir)
				case @settings[:existing_calc_dir]
				when :stop
					Cuby::error("Calculation directory #{calc_dir} already exists")
				when :recreate
					FileUtils.rm_rf(calc_dir)
					calc_dir_create
					Cuby::log.puts_debug("Calculation directory #{calc_dir} recreated")
					return :new
				when :read_results
					Cuby::log.puts_debug("Calculation directory #{calc_dir} kept, results will be reused")
					if output_file
						# Check if output file exists
						unless FileTest.exist?(calc_dir + "/" + output_file)
							Cuby::error("Calculation #{@name} is supposed to reuse existing results\nbut it does not contain output file \"#{output_file}\"") 
						end
					end
					return :old
				when :reuse
					Cuby::log.puts_debug("Calculation directory #{calc_dir} kept, calculation will be run with old input")
					if input_file
						# Check if input file exists
						unless FileTest.exist?(calc_dir + "/" + input_file)
							Cuby::error("Calculation #{@name} is supposed to reuse calculation directory\nbut it does not contain input file \"#{input_file}\"")
						end
					end
					return :old
				end
			else
				Cuby::error("Calculation directory #{calc_dir} can not be created,\nthere is a file with the same name")
			end
		end
		calc_dir_create
		Cuby::log.puts_debug("Calculation directory #{calc_dir} created")
		return :new
	end

	def calc_dir_create
		Dir.mkdir(calc_dir)
		if calc_dir_short?
			# Write long name to the calc dir
			f = File.open(in_calc_dir("_" + @name + "_" + @interface_module::DIRECTORY), "w+")
			f.puts @name + "_" + @interface_module::DIRECTORY
			f.close
		end
	end

	def calc_dir_delete
		#: Delete the calculation directory
		if !FileTest.directory?(calc_dir)
			Cuby::error("Deleting calculation directory:\n\"#{calc_dir}\" does not exist or is not a directory")
		end

		FileUtils.rm_rf(calc_dir)
		Cuby::log.puts_debug("Calculation directory #{calc_dir} deleted")
	end

	CHECKSUM_FILE = '_calc_setup'

	def calc_writing_input_in_calculate
		#: Should be called whenever input is written in calculate, not just in prepare
		#: Prints warning when reused input is overwritten
		if @settings[:existing_calc_dir] == :reuse && @input_reused
			Cuby::warning("New input for calculation \"@name\" must be written but\nkeyword \"existing_calc_dir\" is set to reuse previously prepared one")
		end
	end

	def calc_writing_input
		# Debug message
		Cuby::log.puts_debug("Writing input for calculation #{@name}")

		File.open(calc_dir + "/" + CHECKSUM_FILE, "w+") {|f|
			# Save checksum of current geometry into the calculation directory
			f.puts "geometry_checksum: #{@geometry.checksum}"

			# save @what is calculated
			f.puts "[" + @what.map{|s| ":"+s.to_s}.join(",") + "]"
		}
	end

	def calc_using_input
		# Debug message
		Cuby::log.puts_debug("Calculation #{@name} will use existing input")

		# Set flag on input reuse
		@input_reused = true

		setup_lines = IO.readlines(calc_dir + "/" + CHECKSUM_FILE)

		# Compare checksum of current geometry with the one in directory to be reused,
		# raise error if they do not match
		unless setup_lines[0].strip.split[1] == @geometry.checksum
			Cuby::error("Calculation #{@name}:\nPreviously prepared calculation can not be reused, the geometry is different\n(delete the calculation directories and restart the calculation)")
		end

		# Check @what is calculated
		unless eval(setup_lines[1]) == @what
			Cuby::error("Calculation #{@name}:\nPreviously prepared calculation can not be reused, different properties\nrequested (delete the calculation directories and restart the calculation)")
		end
	end

end
