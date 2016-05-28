
module GeometrySmiles

	def GeometrySmiles.build(smiles_string, settings)

		# Check for presence of the Ballon program
		if settings[:balloon_dir] == "" || !FileTest.exists?(settings[:balloon_dir])
			Cuby::error("To build geometry from the SMILES notation, the Balloon program\n(http://users.abo.fi/mivainio/balloon) must be installed\nand pointed to by the keyword 'balloon_dir'.")
		end

		# Create a temporary filename 
		while FileTest.exist?(filename = 'temp' + sprintf("%08d",(rand*1.0e8).to_i) + ".sdf") do end

		# Run balloon
		balloon_dir = settings[:balloon_dir]
		system "#{balloon_dir}/balloon -f #{balloon_dir}/MMFF94.mff --noGA '#{smiles_string}' #{filename} 2> /dev/null > /dev/null"

		# Test result
		Cuby::error "SMILES translation by Balloon failed:\n#{smiles_string}" unless FileTest.exist?(filename)

		# Read geometry and delete the file
		g = Geometry.new
		g.read_file(filename)
		File.delete(filename)

		return g
	end

end
