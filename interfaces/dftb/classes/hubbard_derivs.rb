module HubbardDerivs

	def HubbardDerivs.from_settings(settings, datafile)
		if settings[:dftb_hubbard_derivs] == :custom
			hubbard_derivs = settings.elements_hash(:dftb_hubbard_derivs_custom)
		else
			File.open(datafile) {|f|
				hubbard_derivs = YAML::load(f)[settings[:dftb_hubbard_derivs]]

			}
		end
		return hubbard_derivs
	end

end

