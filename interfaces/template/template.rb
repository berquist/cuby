################################################################################
#
# Template interface
#
################################################################################

#===============================================================================
# Simplest interface to an external program intended as a template for 
# development of new interfaces.
#
# The calculation setup is simple but the complete logic allowing reusing
# previously prepared calculation directories or just reading of old results
# is fully implemented.
#
# The keywords for this interface are defined in the file 'keywords.yaml'
# located in the same directory
#
# A brief description of the interface provided in 'documentation.html'
# is used in the generation of the manual.
#
# Note on the interface name - it appears at three places:
# 1) name of the directory (under cuby4/interfaces) = inteface name
# 2) name of the main file = interface name + ".rb"
# 3) name of the module = "Interface" + inteface name with first letter capitalized
#===============================================================================

# To read the output of the calculation easily, we will use the parser provided
# by the framework
require "classes/tools/output_parser.rb"

module InterfaceTemplate
	# Each interface is a module that is loaded as an extension of an instance
	# of the Calculation class. It therefore uses some common code defined 
	# in this class.

	#=======================================================================
	# Interface header
	#=======================================================================
	
	# The header defines the properties and capabilities of the interface

	# Development status
	DEVELOPMENT_FLAG = :ok
	DEVELOPMENT_STATUS = "Simplest possible template for interfaces"

	# Interface information
	INTERFACE = :calculation_external # interface using external code for calculation
	CAPABILITIES = [:energy] # It provides only single-point calculation of energy
	MODIFIER = false # It can not be used as a Modifier

	# Suffix used for the calculation directory
	DIRECTORY = "MopacT"

	# Methods provided by the interface:
	# The method must be valid value of the 'method' keyword
	# For each method, available calculated properties are listed separately
	METHODS = {
		:"am1" 		=> [:energy],
		:"pm6" 		=> [:energy],
	}

	#=======================================================================
	# Interface methods
	#=======================================================================
	
	# The following three methods form the API of this interface, these are
	# the entry points called by the framework.

	def prepare_interface
		# This method prepares the calculation. It is called only once
		# at the initializtion of the calculation which can be then
		# called multiple times on changing geometries.

		#...............................................................
		# Prepare calculation directory
		#...............................................................
		
		# If there is existing calculation directory that already
		# contains prepared calculation (the header of the input file)
		# and output of previous calculation, the directory can be
		# reused (depending on the value of 'existing_calc_dir'
		# keyword). Otherwise, new directory is created.

		if calc_dir_mkdir("mopac.header", "mopac.in.out") == :old
			# When the old input is used, the following method
			# is called to check if it is applicable to the current
			# calculation.
			calc_using_input # This method is defined in the Calculation class
			return
		end

		#...............................................................
		# Construction of the input
		#...............................................................

		# It is split into two methods:
	
		# The first one bulids a header file containing the calculation
		# setup (it is called only once)
		mopac_write_keywords # This writes file mopac.header

		# The second one combines the header file with current geometry,
		# producing the final input file.
		mopac_write_input_geo # This writes mopac.in

		#...............................................................
		# Calculation metadata
		#...............................................................

		# Save info on the system upon writing the input - it is checked
		# when existing calculation directory is released.
		calc_writing_input # This method is defined in the Calculation class

		return nil
	end

	def calculate_interface
		# This method executes the calculation on the current geometry
		# It returns the Results object
		
		# The calculation is run unless reading existing results is
		# requested in the input
		mopac_run unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/mopac.in.out")

		# The mopac_read method reads the output file and builds the
		# Results object
		results = mopac_read
		return results
	end

	def cleanup_interface
		# This method cleans up the calulation, in this case just
		# deletes the calculation directory using a shared method
		# provided by the Calculation class.
		calc_dir_delete if @settings[:job_cleanup]
	end

	#=======================================================================
	# Private methods
	#=======================================================================
	
	# These methods do the actual work. They are called from the interface
	# methods defined above but not from outside of this file.

	def mopac_write_keywords
		# Write the header of the input where the calculation is set up
		# geometry is added later, so that the header can be modified
		# by the user after the job is prepared

		# Open the file, the in_calc_dir method prepends the path to the
		# calculation directory
		filename = in_calc_dir("mopac.header")
		f = File.open(filename,"w+")

		# Write the simplest header - method and charge, requesting
		# a single-point calculation
		f.print "#{@settings[:method].upcase} CHARGE=#{@settings[:charge]} 1SCF"

		# Close the file
		f.close
	end

	def mopac_write_input_geo
		# Write final input for mopac, using the previously prepared header
		# and a current geometry

		# Read mopac.header - is just a single line in this case
		header = IO.readlines(in_calc_dir("mopac.header"))[0].strip

		# Open th input file
		filename = in_calc_dir("mopac.in")
		f = File.open(filename,"w+")

		# Print header from the header file, adding option to restart from previous density
		# if it is available and requested in the input.
		if FileTest.exist?(in_calc_dir("mopac.in.den")) && @settings[:start_from_previous]
			restart =  " OLDENS" 
		else
			restart = ""
		end
		f.puts header + restart

		# Separator
		f.print "\n\n"

		# Write the geometry
		# It is in simple element, x, y, z format,
		# default units are Angstroms so no conversion is needed
		@geometry.each_index{|i|
			atom = @geometry[i]
			f.printf("  %2s    %15.10f %15.10f %15.10f\n",atom.element,atom.x,atom.y,atom.z)
		}
		f.puts

		f.close
	end

	def mopac_run
		# Set up environment for MOPAC and run it on the current geometry

		# Update the geometry - it may have changed since the calculation
		# was prepared
		mopac_write_input_geo

		# Prepare the command for running mopac
		# We borrow the keywords fo setup of the external program from the Mopac interface

		# First, go to the calculation directory
		command = "cd #{calc_dir};"
		# Set up the environment - license file is located in the same drectory as the executable
		command << "export MOPAC_LICENSE=#{File.dirname(@settings[:template_mopac_exe])};" 
		# Also, some libraries needed for the run may be located there
		command << "export LD_LIBRARY_PATH=#{File.dirname(@settings[:template_mopac_exe])}:$LD_LIBRARY_PATH;"
		# And finally, the command that runs mopac on the input file
		command << "#{@settings[:template_mopac_exe]} mopac.in 2> mopac.err;"

		# Run mopac, exit if it returned bad exit code
		unless system(command)
			Cuby::error "Mopac returned nonzero exit code (calculation #{@name})"
		end
	end

	def mopac_read 
		# Read the output produced by successful MOPAC run
		# Returns new instance of the Results class
	
		# Create empty Results object
		results = Results.new

		# Create parser acting on the ouput file
		parser = OutputParser.new(in_calc_dir("mopac.in.out"), @name)

		# Parser pattern for energy
		# The arguments are:
		# 1) name of the result
		# 2) regular expression to look for
		# 3) further options how to handle the line - in this case, the submatch, it is the part of
		#    the regexp in brackets is saved - this is the actual number
		parser.add_pattern(:energy, /^ +FINAL HEAT OF FORMATION = +([-]*[0-9]*\.[0-9]+)/, :get => :submatch)

		# Errors
		# The parser can also look for error messages in the output and translate them to Cuby
		# errors
		parser.add_error(/CHARGE SPECIFIED FOR SYSTEM .* IS INCORRECT/, "Mopac thinks the charge is wrong")


		# Run parser
		parser.execute

		# Save the value extracted by the parser to the Results object
		# The result is in kcal/mol, no conversion is needed
		results.energy = parser[:energy]

		return results
	end

end
