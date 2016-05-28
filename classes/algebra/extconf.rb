require 'mkmf'
require 'yaml'

# Configure paths to external libraries
# (optional, if the libraries are not at a place recognized by the system)
fn = File.expand_path("~/.cuby4")
if File.exists?(fn)
	File.open(fn){|f|
		yamldata = YAML.load(f)
		if @libdirs = yamldata["extension_lib_dirs"]
			@libdirs.each{|dir|
				dir_config('',dir)
			}
		end
	}
end

# Search for library file in libdirs
def lib_file_found(libname)
	@libdirs.each{|dir|
		Dir.foreach(dir){|file|
			if file == "lib"+libname+".so"
				puts "Library #{libname} found, as a file in #{dir}"
				return true
			end
		}
	}
	return false
end


# Blas/Lapack: needed for all computationally-intensive matrix operations

unless ARGV.include?("--no-blas")
	if ARGV.include?("--force-blas")
		$CFLAGS += " -DBLAS_FOUND" # Define flag used in C code
		$LIBS += " -lblas"
	elsif  have_library("blas")
		$CFLAGS += " -DBLAS_FOUND" # Define flag used in C code
	elsif lib_file_found('blas')
		$CFLAGS += " -DBLAS_FOUND" # Define flag used in C code
		$LIBS += " -lblas"
	else
		puts "BLAS library not found, extension will be compiled without BLAS support"
	end
end

unless ARGV.include?("--no-lapack")
	if ARGV.include?("--force-lapack")
		$CFLAGS += " -DLAPACK_FOUND" # Define flag used in C code
		$LIBS += " -llapack" # Add directly to libs
	elsif have_library("lapack")
		$CFLAGS += " -DLAPACK_FOUND" # Define flag used in C code
	elsif lib_file_found('lapack')
		$CFLAGS += " -DLAPACK_FOUND" # Define flag used in C code
		$LIBS += " -llapack" # Add directly to libs
	else
		puts "LAPACK library not found, extension will be compiled without some routines"
	end
end

# UMFPACK: sparse matrix solver - optional
unless ARGV.include?("--no-umfpack")
	if have_library("umfpack")
		puts "AMD library not found" unless have_library("amd")
		$CFLAGS += " -DUMFPACK_FOUND" # Define flag used in C code
	else
		puts "UMFPACK library not found, extension will be compiled without sparse matrix solver"
	end
end

# Optimization
$CFLAGS += " -O3"

# Create makefile
create_makefile("algebra_c")

# Save ruby version used for compilation
File.open("_version","w+"){|f|
	f.puts RUBY_VERSION
}
