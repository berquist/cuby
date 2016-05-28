require "classes/objects/gradient.rb"
require "classes/objects/hessian.rb"
require "classes/objects/atomic_charges.rb"
require "classes/objects/multipoles.rb"
require "classes/objects/polarizability.rb"
require "classes/objects/molecular_orbitals.rb"

class Results
	attr_accessor :energy
	attr_accessor :energy_components

	attr_accessor :gradient
	attr_accessor :gradient_components

	attr_accessor :hessian

	attr_accessor :atomic_charges
	attr_accessor :atomic_charges_gradient

	attr_accessor :multipoles
	attr_accessor :static_polarizability

	attr_accessor :point_charges_gradient

	attr_accessor :molecular_orbitals

	# energy_type :total, :interaction
	# optimizer status
	# optimizer_cycles
	# md_cycles
	# md_restart_data (velocities ...)

	def initialize
		@energy = nil
		@energy_components = {}

		@gradient = nil
		@gradient_components = {}

		@hessian = nil

		@atomic_charges = nil

		@multipoles = Multipoles.new # Empty container
		@static_polarizability = nil
	end

	def add_modifier_results(modifier)
		@energy += modifier.results.energy
		@energy_components[('modifier_' + modifier.settings[:interface].to_s).to_sym] = modifier.results.energy

		if @gradient
			@gradient.plus!(modifier.results.gradient)
		end

		if @hessian
			@hessian += modifier.results.hessian
		end
	end

	#=======================================================================
	# Printing
	#=======================================================================

	def print_energy(name = "Energy")
		Cuby::log.puts_v_only(:minimal, '%.6f' % @energy)
		Cuby::log.puts_v_only([:brief, :normal], "#{name}: #{'%.6f' % @energy} kcal/mol")
	end

	def print_energy_decomposition(name = "Energy")
		Cuby::log.puts_v_only([:brief, :normal], name + " decomposition:")
		Cuby::log.puts_v_only([:brief, :normal], "   no components found") if @energy_components.empty?
		@energy_components.each_pair{|key, value|
			Cuby::log.puts_v_only(:minimal, sprintf("%s %.6f",key.to_s, value))
			Cuby::log.puts_v_only([:brief, :normal], sprintf("   %-30s%16.6f kcal/mol",key.to_s, value))
		}
	end

	def print_atomic_charges(geometry)
		Cuby::log.puts_v_only([:brief, :normal], "Atomic charges (#{@atomic_charges.type_s}):")
		geometry.each_index{|i|
			Cuby::log.puts(sprintf("%4s%10.4f",geometry[i].element, @atomic_charges[i]))
		}
	end

	def print_dipole
		return unless @multipoles
		return unless @multipoles[:dipole]
		Cuby::log.puts @multipoles[:dipole]
	end

	def print_gradient
		Cuby::log.puts_v_only([:brief, :normal], "Gradient (kcal/mol/A):")
		Cuby::log.puts_v(:minimal, @gradient.to_s)
	end

	def print_polarizability
		Cuby::log.puts_v(:minimal, @static_polarizability) if @static_polarizability
	end

	def print_mos
		return unless @molecular_orbitals
		Cuby::log.puts_v_only([:brief, :normal], "\nMolecular orbitals:")
		@molecular_orbitals.each{|mo|
			Cuby::log.puts_v(:minimal, mo)
		}
		Cuby::log.puts_v(:minimal, "HOMO energy #{@molecular_orbitals.homo.energy_au} a.u.")
		Cuby::log.puts_v(:minimal, "LUMO energy #{@molecular_orbitals.lumo.energy_au} a.u.")
	end
end
