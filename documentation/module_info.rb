
module ModuleInfo
	def self.all_keywords(module_dir)
		# Search module directory for all keywords used

		# List of files in the module directory and its subdirectories
		# (one level only)
		filelist = []
		Dir.entries(module_dir).sort.each{|dir|
			next if dir == ".."
			dir = module_dir + "/" + dir
			if FileTest.directory?(dir)
				Dir.entries(dir).sort.each{|fn|
					next if fn == ".." || fn == "."
					fn = dir + "/" + fn
					next if FileTest.directory?(fn)
					filelist << fn.gsub("/./","/")
				}
			end
		}

		# Add cuby files loaded by these files
		filelist2 = []
		filelist.each{|fn|
			# Ruby files only
			next unless fn =~ /\.rb$/
			filelist2 |= ModuleInfo.scan_required_cuby_files(fn)
		}
		filelist |= filelist2

		# Scan the files for keywords
		keywords = []
		filelist.each{|fn|
			# Ruby files only, including ERB templates
			next unless fn =~ /\.rb$/ || fn =~ /\.erb$/
			keywords |= ModuleInfo.scan_file_for_keywords(fn)
		}
		return keywords.sort
	end

	def self.scan_required_cuby_files(fn)
		cuby_dir = File.expand_path(File.dirname(__FILE__)+"/..")
		required = []
		IO.readlines(fn).each{|line|
			begin
				next if line =~ /^\s*#/ # skip comments
			rescue
				$stderr.puts "ASCII encoding error in file #{fn}"
				exit
			end
			# Shared code in cuby should be located in classes directory
			if match = /require *['"](classes.*)['"]/.match(line)
				path = cuby_dir + "/" + match[1]
				if FileTest.file?(path)
					required << path
				elsif FileTest.file?(path + ".rb")
					required << path + ".rb"
				end
			end
		}
		return required.uniq
	end

	def self.scan_file_for_keywords(fn)
		keywords = []
		IO.readlines(fn).each{|line|
			begin
				next if line =~ /^\s*#/ # skip comments
			rescue
				$stderr.puts "ASCII encoding error in file #{fn}"
				exit
			end
			line.scan(/settings\[:([^\]]+)\]/){|kw|
				kw2 = kw[0].split(/\s*,\s*:/).last
				if Settings.all_keywords[kw2.to_sym]
					keywords << kw2 unless Settings.all_keywords[kw2.to_sym].devel
				end
			}
			line.scan(/settings.elements_hash\(:([^\]]+)\)/){|kw|
				kw2 = kw[0].split(/\s*,\s*:/).last
				if Settings.all_keywords[kw2.to_sym]
					keywords << kw2 unless Settings.all_keywords[kw2.to_sym].devel
				end
			}
		}
		return keywords.uniq
	end

	def self.local_keywords(module_dir)
		fn = module_dir + "/" + "keywords.yaml"
		return [] unless FileTest.exist?(fn)
		keywords = []
		IO.readlines(fn).each{|line|
			line.scan(/:([^:]+): *!ruby\/object:SettingsKeyword/){|kw|
				keywords << kw[0] unless Settings.all_keywords[kw[0].to_sym].devel
			}
		}
		return keywords.sort
	end
end
