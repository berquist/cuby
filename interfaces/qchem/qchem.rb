################################################################################
#
# Q-Chem interface
#
# Author: Eric Berquist
# Date created: 2016-05-29
# License: Cuby4 license
# Description: Interface for external calculations using Q-Chem
# Status: Works
#
################################################################################

#===============================================================================
# Interface to Q-Chem
# http://www.q-chem.com/
#===============================================================================

require "classes/tools/output_parser.rb"

module InterfaceQchem
    #=======================================================================
    # Interface header
    #=======================================================================
    DEVELOPMENT_FLAG = :ok
    DEVELOPMENT_STATUS = "OK"
    # Interface information
    INTERFACE = :calculation_external
    CAPABILITIES = [:energy]
    MODIFIER = false
    DIRECTORY = "QChem"
    # Methods provided by the interface:
    METHODS = {
        :"hf"       => [:energy]
    }
    # Input structure
    INPUT_BLOCKS = [
        InputBlock[:molecule_a, :optional, "Definition of the first monomer (selection, charge, multiplicity)"],
        InputBlock[:molecule_b, :optional, "Definition of the second monomer (selection, charge, multiplicity)"]
    ]
    #=======================================================================

    def prepare_interface
        # Prepare calculation directory
        if calc_dir_mkdir("qchem.header", "qchem.out") == :old
            calc_using_input # Old input should be reused
            return
        end

        # Write the keyword part of the input
        qchem_write_keywords # This writes file qchem.header

        # Create complete input by combining qchem.header and the current geometry
        qchem_write_input_geo

        # Save info on the system upon writing the input
        calc_writing_input
    end

    def calculate_interface
        qchem_run unless @settings[:existing_calc_dir] == :read_results && FileTest.exist?(calc_dir+"/qchem.out")
        return qchem_read
    end

    def cleanup_interface
        calc_dir_delete if @settings[:job_cleanup]
    end

    #=======================================================================
    # Private methods
    #=======================================================================

    def qchem_write_keywords
        filename = in_calc_dir("qchem.header")
        f = File.open(filename,"w+")

        f.puts "$rem"

        f.puts " method #{@settings[:method]}"
        # when :dft
        #     require "classes/calculation/dft_options.rb"
        #     dft_file = File.open(interface_dir + "/dft_options.yaml", "r")
        #     dft_options = YAML.load(dft_file)
        #     dft_file.close
        #     f.puts "#{what}('#{dft_options.functional(@settings)}')"

        f.puts " basis #{@settings[:basisset]}"
        if @settings.set?(:auxiliary_basis_mp2)
            f.puts " aux_basis #{@settings[:auxiliary_basis_mp2]}"
        end

        if @settings[:spin_restricted]
            f.puts " unrestricted false"
        else
            f.puts " unrestricted true"
        end

        if @settings[:correlation_frozen_core]
            f.puts " n_frozen_core fc"
        else
            f.puts " n_frozen_core 0"
        end

        f.puts " scf_convergence  #{@settings[:scf_convergence]}" if @settings.set?(:scf_convergence)

        f.puts " symmetry false" unless @settings[:use_symmetry]
        f.puts " sym_ignore true" unless @settings[:use_symmetry]
        f.puts " mem_static #{@settings[:mem]}"
        # How to handle `mem_total` and `cc_memory`?

        f.puts "$end"

        f.puts "$molecule"
        f.puts "%geometry%"
        f.puts "$end"

        f.close
    end

    def qchem_write_input_geo
        # Write input file using the existing header, only geometry is added
        filename = in_calc_dir("qchem.in")
        f = File.open(filename,"w+")

        # Template from file, substitute geometry
        IO.readlines(in_calc_dir("qchem.header")).each{|line|
            if line =~ /%geometry%/
                # Write geometry')"

                if @settings[:qchem_geometry_fragments]
                    # Fragments written separately
                    qchem_write_geo_fragments(f)
                else
                    qchem_write_geo_block(f, @geometry, @settings[:charge], @settings[:multiplicity])
                end
            else
                # Write the original line
                f.puts line
            end
        }

        f.close
    end

    def qchem_write_geo_fragments(file)
        # Assumes there are two blocks with molecule definition
        [:molecule_a, :molecule_b].each{|molname|
            file.puts "--" unless molname == :molecule_a
            geo = @geometry.geometry_from_selection(@settings[molname, :selection])
            charge = @settings[molname, :charge]
            multiplicity = @settings[molname, :multiplicity]
            qchem_write_geo_block(file, geo, charge, multiplicity)
        }
    end

    def qchem_write_geo_block(file, geo, charge, multiplicity)
        file.puts "#{charge} #{multiplicity}"
        geo.each_index{|i|
            atom = geo[i]
            element = atom.element.to_s
            element = "@" + element if atom.properties[:ghost]
            file.printf("%5s%25.15f%25.15f%25.15f\n", element, atom.x, atom.y, atom.z)
        }
    end

    def qchem_run
        # Build input with the current geometry
        qchem_write_input_geo

        # Run qchem
        command = "cd #{calc_dir};"
        command << "$(which qchem) -nt #{@settings[:parallel]} qchem.in qchem.out > output.yaml 2> qchem.err;"
        unless system(command)
            # Generic error
            Cuby::error "Q-Chem returned nonzero exit code (calculation #{@name})"
        end
    end

    def qchem_read # => Results
        results = Results.new
        parser = OutputParser.new(in_calc_dir("qchem.out"), @name)

        parser.add_pattern(:e_scf, /^\s*Total energy in the final basis set = /)

        parser.execute

        case @settings[:method]
        when :hf
            results.energy = parser[:e_scf] * HARTREE2KCAL
        end

        return results
    end

end
