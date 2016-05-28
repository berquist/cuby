
module SearchForFile
	def SearchForFile.search(arg_hash)
		# Default arguments
		args = {
			:locations => ['./', ''],
			:extensions => [''],
			:raise_error => true
		}
		args.merge!(arg_hash) # Overwrite defaults as needed

		# Check for filename
		raise "SearchForFile.search requires argument :filename" unless args[:filename]


		# Scan the locations
		args[:locations].each{|loc|
			loc = loc + "/" unless loc =~ /\/$/
			# Scan the extensions
			args[:extensions].each{|ext|
				fn = loc + args[:filename] + ext
				# expand path, ruby does not recognize e.g. "~/"
				fn = File.expand_path(fn)
				return fn if FileTest.exist?(fn)
			}
		}

		if args[:raise_error]
			# Error message if not found
			files = []
			args[:extensions].each{|ext|
				files << args[:filename] + ext
			}
			err_str = "File #{files.join('|')} not found in:"
			args[:locations].each{|loc|
				err_str += "\n   #{loc}"
			}
			Cuby::error(err_str)
		else
			# Return false
			return nil
		end
	end
end

