require "classes/objects/results.rb"
require "classes/objects/point_charges.rb"
require "classes/calculation/calculation_queue.rb"
require "classes/calculation/calculation_remote.rb"
require "classes/calculation/default_interface_selection.rb"
require "classes/calculation/input_block.rb"

class CalculationError < RuntimeError
end

class Calculation
	include CalculationRemote

	# Attributes
	attr_reader :name
	attr_reader :geometry
	attr_reader :settings
	attr_reader :queue
	attr_reader :what

	attr_accessor :point_charges

	# Optionally, protocols may refence themselves here:
	attr_accessor :parent_protocol

	# Allowed interface types
	INTERFACE_TYPES = [
		:calculation,
		:calculation_external,
		:composite
	]

	@@interfaces = {} # List of available interfaces
	@@names = [] # List of calculation names (To ensure they are unique)

	#=======================================================================
	# Build list of interfaces available in the interfaces directory
	#=======================================================================
	@@interface_path = Cuby.install_dir + "/interfaces"
	# Fill the list
	Dir.entries(@@interface_path).sort.each{|dir|
		next if dir == '.' || dir == '..'
		next unless FileTest.directory?(@@interface_path + "/" + dir)
		@@interfaces[dir.downcase.to_sym] = @@interface_path + "/" + dir
	}

	#=======================================================================
	# Access to the list of interfaces and interface information
	#=======================================================================
	def Calculation.available_interfaces
		return @@interfaces
	end

	#=======================================================================
	# Constructor
	#=======================================================================
	
	def initialize(name, settings, geometry)
		# Check for unique name
		if @@names.include?(name)
			raise(CalculationError, "Each instance of Calculation must have unique name (\"#{name}\" already exists)")
		end
		@@names << name

		@name = name
		@geometry = geometry
		@settings = settings
		@results = nil
		@point_charges = nil
		@modifiers = []

		# Prepare point charges from settings
		if @settings[:point_charges] != ""
			@point_charges = PointCharges.load_from_settings(@settings[:point_charges])
			Cuby::log.puts_debug("Loading point charges from settings (calculation #{name}, #{@point_charges.size} charges)")
		end
		
	end

	#=======================================================================
	# Helpers
	#=======================================================================
	
	def interface_dir
		return @@interface_path + "/" + @settings[:interface].to_s
	end
	
	#=======================================================================
	# Methods
	#=======================================================================

	def prepare(what, option = nil)
		@what = what

		# Interface from method
		if @settings[:interface] == :auto
			@settings[:interface] = DefaultInterfaceSelection.from_settings(@settings)
		end

		# Check whether interface exists
		unless @@interfaces[@settings[:interface]]
			Cuby::error("Requested interface \"#{@settings[:interface]}\" does not exist")
		end

		# Check if there are keywords for other interfaces
		@settings.each_keyword{|kw|
			kwd = Settings::keyword(kw)
			if kwd.interface && ! kwd.interface.include?(@settings[:interface])
				s = kwd.interface.size > 1 ? "s" : ""
				Cuby::warning "Keyword '#{kw}' ignored in calculation '#{@name}', it belongs to other interface#{s} (#{kwd.interface.join(', ')})"
			end
		}

		# Load interface
		interface_dir =  @@interfaces[@settings[:interface]]
		interface_fn = interface_dir + "/" + @settings[:interface].to_s + ".rb"
		unless FileTest.exist?(interface_fn)
			Cuby::error("Interface \"#{@settings[:interface]}\": file \"#{@settings[:interface].to_s}.rb\" not found in the interface directory")
		end
		require interface_fn
		@interface_fn = interface_fn
		@interface_module = eval("Interface" + @settings[:interface].to_s.split("_").map{|s| s.capitalize}.join(''))

		# Interface check: interface type
		unless INTERFACE_TYPES.include?(@interface_module::INTERFACE)
			Cuby::error("Interface \"#{@settings[:interface]}\" has wrong type (defined in INTERFACE constant)")
		end
		# Interface check: development status (raise warning when interface in development is called)
		case @interface_module::DEVELOPMENT_FLAG
		when :ok
			# Nothing to be done
		when :warning
			Cuby::warning("Interface #{@settings[:interface]} is under development and should be used with care:\n"+@interface_module::DEVELOPMENT_STATUS, :once)
		when :error
			Cuby::error("Interface #{@settings[:interface]} is under development and can not be used")
		else
			raise "Interface #{@settings[:interface]} does not have proper DEVELOPMENT_FLAG constant"
		end

		# Capabilities check: is the interface a modifier?
		if option == :modifier
			unless @interface_module::MODIFIER
				Cuby::error("Interface \"#{@settings[:interface]}\" can not be used as a modifier")
			end
		end

		# Capabilities check: properties inplemented
		ignored_tasks = []
		@what.each{|task|
			if @interface_module::CAPABILITIES.include?(task)
				# Everything OK
			elsif @interface_module.const_defined?('IGNORE') && @interface_module::IGNORE.include?(task)
				# Some options are ignored
				ignored_tasks << task
			elsif option == :modifier && @interface_module.const_defined?('IGNORE_AS_MODIFIER') && @interface_module::IGNORE_AS_MODIFIER.include?(task)
				# When run as modifier, some tasks are ignored
				ignored_tasks << task
			else
				# Not available, raise error
				Cuby::error("Interface \"#{@settings[:interface]}\" does not support calculation of #{task}")
			end
		}

		# Capabilities check: ghost and dummy atoms
		if !@interface_module::CAPABILITIES.include?(:ghost_atoms) && @geometry.ghost_atoms? && !ignored_tasks.include?(:ghost_atoms)
			Cuby::error("Interface \"#{@settings[:interface]}\" does not support ghost atoms")
		end
		if !@interface_module::CAPABILITIES.include?(:dummy_atoms) && @geometry.dummy_atoms? && !ignored_tasks.include?(:dummy_atoms)
			Cuby::error("Interface \"#{@settings[:interface]}\" does not support dummy atoms")
		end
		# Capabilities check: electric field
		if @settings.set?(:electric_field) && !@interface_module::CAPABILITIES.include?(:electric_field) && !ignored_tasks.include?(:electric_field)
			Cuby::error("Interface \"#{@settings[:interface]}\" does not support external electric field")
		end
		# Capabilities check: point_charges
		if @point_charges && !@interface_module::CAPABILITIES.include?(:point_charges) && !ignored_tasks.include?(:point_charges)
			Cuby::error("Interface \"#{@settings[:interface]}\" does not support external point charges")
		end

		# Capabilities check: atomic charges
		if @what.include?(:atomic_charges) && !ignored_tasks.include?(:atomic_charges)
			begin
			       	@interface_module::ATOMIC_CHARGES
			rescue
				Cuby::error("Interface \"#{@settings[:interface]}\" does not support atomic charges calculations")
			end
			unless @interface_module::ATOMIC_CHARGES.include?(@settings[:atomic_charges])
				Cuby::error("Interface \"#{@settings[:interface]}\" does not support atomic charges model \"#{@settings[:atomic_charges]}\"")
			end
		end

		# Capabilities check: Solvent
		if @settings[:solvent_model] != :none && !ignored_tasks.include?(:solvent)
			begin
			       	@interface_module::SOLVENTS
			rescue
				Cuby::error("Interface \"#{@settings[:interface]}\" does not support calculations in solvent")
			end
			unless @interface_module::SOLVENTS.include?(@settings[:solvent_model])
				Cuby::error("Interface \"#{@settings[:interface]}\" does not support solvent model \"#{@settings[:solvent_model]}\"")
			end
		end


		Cuby::log.puts_debug("Preparing calculation \"#{@name}\"")
		Cuby::log.puts_debug("   interface: #{@settings[:interface]}")
		Cuby::log.puts_debug("   calculate: #{@what.map{|x| x.to_s}.join(', ')}")

		# Extend this object with interface module
		extend @interface_module

		# Check mandatory methods
		if @interface_module::INTERFACE == :composite
			method_list = [:prepare_interface, :queue_interface, :compose_interface]
		else
			method_list = [:prepare_interface, :calculate_interface]
		end
		method_list.each {|method|
			unless self.respond_to?(method)
				Cuby::error("Interface \"#{@settings[:interface]}\" does not have method #{method}")
			end
		}

		if @interface_module::INTERFACE == :calculation_external
			# Load methods shared by external calculations
			require "classes/calculation/calculation_external.rb"
			extend CalculationExternal

			# Some interfaces provide multiple methods: check whether the requested method is available
			unless @settings.set?(:method)
				if @interface_module::METHODS.size == 1
					# The single method is used by default
					@settings[:method] = @interface_module::METHODS.keys[0]
					Cuby::log.puts_debug "Default method '#{@settings[:method]}' used. (#{@name})"
				else
					Cuby::error("The interface provides more methods, select one using the 'method' keyword.\n(#{@name})")
				end
			end
			# Check whether the interface provides the requested method
			check_method_availability
			# Check method_capabilities - e.g. gradient may not be available for all methods
			check_method_capabilities
		end


		# Prepare interface
		prepare_interface

		# Prepare modifiers
		# Build calculations
		if @settings.set?(:modifiers)
			@settings[:modifiers].each{|mod_name|
				# Check whether interface exists
				unless @@interfaces[mod_name]
					Cuby::error("Interface \"#{mod_name}\" requested as a modifier does not exist")
				end

				# Look for settings block
				blockname = ("modifier_" + mod_name.to_s).to_sym
				unless @settings.block(blockname)
					# Create empty if none found
					 @settings.new_block(blockname)
				end
				mod_settings = @settings.block(blockname)
				# Check for interface match
				if mod_settings.set?(:interface)
					if mod_settings[:interface] != mod_name
						Cuby::error("Modifier interface (\"#{mod_settings[:interface]}\") does not match modifier name (\"#{mod_name}\")")
					end
				end

				# Set interface
				mod_settings[:interface] = mod_name
				# Set required settings

				@modifiers << Calculation.new(@name + '_mod_' + mod_name.to_s, mod_settings, @geometry)
			}

			# Run prepare
			@modifiers.each{|modifier|
				modifier.prepare(@what, :modifier)
			}
		end
	end

	def send_to_queue(queue_obj)
		# Queue this
		@results = nil
		@results_processed = false
		@queue = queue_obj

		if @interface_module::INTERFACE == :composite
			# Interface may have subcalculations that are queued separately
			queue_interface(queue_obj)
		else
			# Otherwise, queue this
			queue_obj.add_to_queue(self)
		end

		# Call queue on modifiers
		@modifiers.each{|calc| calc.send_to_queue(queue_obj)}
	end

	def calculate(thread_no)
		# Called by the queue system
		time = Time.now
		if @interface_module::INTERFACE == :calculation || @interface_module::INTERFACE == :calculation_external
			case @@remote
			when :no
				# Call the interface directly
				@results = calculate_interface
			when :drb
				# Call it on the remote server
				calc_server = @@remote_drb_servers[thread_no]
				wd = remote_copy_workdir_to_server(calc_server.hostname, calc_server.workdir) if @interface_module::INTERFACE == :calculation_external
				@results = calc_server.calculate(Marshal.dump(self), @interface_fn, wd) # For some reason, the calculation object must be serialized
				remote_copy_workdir_from_server(calc_server.hostname, calc_server.workdir) if @interface_module::INTERFACE == :calculation_external
			end
		else
			raise "Method 'calculate' can be called only for interfaces of type ::calculation or :calculation_external"
		end
		Cuby::log.puts_debug("Calculation \"#{@name}\" took #{'%.2f' % (Time.now - time)} seconds")
	end

	def results
		# There are no results when only prepare was run, do not try to get them
		return nil if @settings[:prepare_only]

		# If no results exist
		unless @results
			# Trigger queue execution
			@queue.execute 
		end

		unless @results_processed
			@results_processed = true

			# Composite interfces may need to process the components and build the final results
			if @interface_module::INTERFACE == :composite
				@results = compose_interface
			end

			# Apply_modifiers
			@modifiers.each{|modifier|
				@results.add_modifier_results(modifier)
			}

			# Save lookup
			if @settings[:calculation_print_lookup]
				puts "LOOKUP>  #{@geometry.checksum}:" 
				puts "LOOKUP>    energy: #{@results.energy}"
			end

			# Replace the energy with a component if desired
			if @what.include?(:energy) 
				if @settings.set?(:replace_energy_with_component) && @settings.set?(:replace_energy_with_expression)
					# Warning when both options are used
					Cuby::warning "Both replace_energy_with_expression and replace_energy_with_component are set,\n replace_energy_with_expression is used."
				end

				if @settings.set?(:replace_energy_with_expression)
					@results.energy_components[:original_energy] = @results.energy
					expression = " " + @settings[:replace_energy_with_expression] + " "
					# Replace names of components in the expression with ruby variables
					@results.energy_components.each_pair{|key, value|
						expression = expression.gsub(
							/([^a-z,A-Z,0-9,_])#{key.to_s}([^a-z,A-Z,0-9,_])/,
							"\\1@results.energy_components[:#{key.to_s}]\\2"
						)
					}
					@results.energy = eval(expression)
				elsif @settings.set?(:replace_energy_with_component)
					@results.energy_components[:original_energy] = @results.energy
					component = @settings[:replace_energy_with_component]
					Cuby::error "Component '#{component}' that should be used instead of energy not found" unless @results.energy_components[component]
					@results.energy = @results.energy_components[component]
				end
			end
		end


		# return results
		return @results
	end

	def cleanup
		# Clean this
		cleanup_interface if self.respond_to?(:cleanup_interface)
		# Clean modifiers
		@modifiers.each{|modifier|
			modifier.cleanup
		}
		# Remove from list of names
		@@names.delete(@name)
	end

	def priority
		return priority_interface if self.respond_to?(:priority_interface)
		return @geometry.size.to_f
	end
end

