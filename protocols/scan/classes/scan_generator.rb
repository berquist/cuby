class ScanGenerator
	def check_setup(settings, keys)
		keys.each{|key|
			unless settings[:scan_generator_setup].has_key?(key)
				Cuby::error "Missing key '#{key}' in keyword 'scan_generator_setup'"
			end
		}
	end

	def range_from_value(value, keyname)
		if value.kind_of? Numeric
			return [value.to_f]
		elsif value.kind_of? Array
			return value
		elsif value.kind_of? String
			return Settings::Seq.from_string(value).to_array
		else
			Cuby::error "Key '#{key}' in scan_generator_setup has to be a number, array or a sequence command"
		end
	end

	def run(geometry, output_fn)
		f = File.open(output_fn, "w+")
		count = iterate(geometry){|geo, str|
			geo.write_xyz(:file => f, :second_line => str)
		}
		f.close
		Cuby::log.puts "Generated #{count} points"
		return nil
	end

	def build_batches(geometry, batchsize, dirname_prefix, filename)
		i = 0
		f = nil
		digits = 5
		dirlist = []
		count = iterate(geometry){|geo, str|
			if i % batchsize == 0
				f.close if f
				dirname = sprintf("#{dirname_prefix}_%0#{digits}d", i)
				Dir.mkdir(dirname)
				dirlist << dirname
				f = File.open(dirname+"/"+filename, "w+")

			end
			geo.write_xyz(:file => f, :second_line => str)
			i += 1
		}
		f.close
		Cuby::log.puts "Generated #{count} points"
		Cuby::log.puts "Split into #{dirlist.size} batches"
		Cuby::log.puts ""
		return dirlist
	end
end
