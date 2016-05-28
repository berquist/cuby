require "classes/lib/minitar"
require "stringio"

module GeometryDatabase
	include Archive::Tar

	DB_PATH = "/data/geometries"
	DB_PATH_DEV = "/data/geometries_devel"

	# Cache for database tables
	@@databases = {}	# Database data
	@@databases_paths = {}	# Database locations

	def GeometryDatabase.get(name, settings)
		# Returns a geometry loaded from the database.
		# The item is specified by a name in following format:
		# "db_name:item_name"

		# Check syntax
		unless name =~ /[^:]+:[^:]/
			Cuby::error "Wrong format of the geometry database name 'name'"
		end

		# Get DB name and item name
		db, item = name.strip.split(":")
		item.downcase!

		item = item.to_i if item =~ /^[0-9]+$/

		# Look for the DB
		if @@databases.has_key?(db)
			# database is already cached
		elsif File.exists?(Cuby::install_dir + DB_PATH + "/#{db}.yaml")
			# database file exists in default location - load into cache now
			GeometryDatabase.load_file(db, DB_PATH)
		elsif  File.exists?(Cuby::install_dir + DB_PATH_DEV + "/#{db}.yaml")
			# database file exists in development location - load into cache now
			GeometryDatabase.load_file(db, DB_PATH_DEV)
			Cuby::warning("Geometry database loaded from a development location")
		elsif File.exists?(settings[:geometry_database_path] + "/#{db}.yaml")
			# database file exists at user location - load into cache now
			GeometryDatabase.load_file(db, settings[:geometry_database_path])
		else
			Cuby::error "Geometry database '#{db}' does not exist"
		end

		# Look for the geometry
		unless filename = @@databases[db][item]
			Cuby::error "Geometry '#{item}' not found in database '#{db}'"
		end

		# Load the geometry file from tar archive
		geo = GeometryDatabase.read_tar(@@databases_paths[db] + '/' + db + '.tar', filename)

		# Return geometry
		return geo
	end

	def GeometryDatabase.read_tar(archivename, filename)
		# Open the archive file
		begin
			a = File.open(archivename, "rb")
		rescue
			Cuby::error "Error opening geometry database archive '#{archivename}'"
		end

		begin
			# Search for the entry
			Minitar::Input.open(a) {|inp|
				inp.each {|entry|
					if entry.full_name == filename
						# Geometry file found, read it
						filetype = GeoFile.type_from_filename(filename)
						f = StringIO.new(entry.read)
						geo = Geometry.new
						GeoFile::read(geo, filetype, :file => f)
						f.close
						return geo
					end
				}
			}
			a.close
		rescue
			Cuby::error "Error reading geometry '#{filename}' from database archive '#{archivename}'"
		end

		# Geometry not found, raise error
		Cuby::error "Geometry '#{filename}' not found in database archive '#{archivename}'"
	end

	#=======================================================================
	# Private
	#=======================================================================
	
	def GeometryDatabase.load_file(db, path)
		@@databases_paths[db] = Cuby::install_dir + path
		filename = Cuby::install_dir + path + "/" + db + ".yaml"
		File.open(filename) {|f|
			@@databases[db] = {}
			YAML::load(f).each_pair{|key,val| 
				key = key.downcase if key.kind_of?(String)
				key = key.to_i if key =~ /^[0-9]+$/
				@@databases[db][key] = val
			}
		}
	end

end
