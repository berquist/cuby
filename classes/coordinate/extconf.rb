require 'mkmf'

# Optimization
$CFLAGS += " -O3"

create_makefile("coordinate_c")

# Save ruby version used for compilation
File.open("_version","w+"){|f|
	f.puts RUBY_VERSION
}
