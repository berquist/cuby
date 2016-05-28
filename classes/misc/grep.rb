
require 'iconv' unless String.method_defined?(:encode) # For unicode-safe version

module Grep
	def Grep.file(filename, pattern)
		File.open(filename) {|f|
			while line = f.gets
				return true if line =~ pattern
			end
		}
		return false
	end

	def Grep.file_unicode(filename, pattern)
		unless String.method_defined?(:encode)
			ic = Iconv.new('UTF-8', 'UTF-8//IGNORE')
		end
		File.open(filename) {|f|
			while line = f.gets
				if String.method_defined?(:encode)
					line.encode!('UTF-8', 'UTF-8', :invalid => :replace)
				else
					line = ic.iconv(line)
				end
				begin
					return true if line =~ pattern
				rescue
					# Rescue any error caused by invalid characters
				end
			end
		}
		return false
	end
end
