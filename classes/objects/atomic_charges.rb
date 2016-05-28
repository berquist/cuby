
# Array of atomic charges, indexes correspond to indexes in a Geometry

class AtomicCharges < Array
	CHARGE_TYPES = {
		:unknown	=> "Charges from unknown source",
		:mulliken	=> "Mulliken charges",
		:eem		=> "Electrostatic Equilibration Method (empirical) charges",
		:spatial	=> "Spatial population charges",
		:loewdin	=> "Loewdin charges",
		:becke		=> "Becke charges",
		:hirshfeld	=> "Hirshfeld charges",
		:nbo		=> "Natural bond orbital (NBO) charges",
		:forcefield	=> "Forcefield charges"
	}

	attr_reader :charge_type
	attr_reader :gradient

	def initialize(type, size = 0)
		check_charges_type(type)
		@charge_type = type
		super(size, 0.0)
	end

	def check_charges_type(type)
		unless CHARGE_TYPES.has_key?(type)
			raise "Unknown type of atomic charges: '#{type}'"
		end
	end

	def type_s
		return CHARGE_TYPES[@charge_type]
	end

	def write_to_geometry(geometry)
		unless geometry.size == self.size
			raise "Number of atomic charges does not match number of atoms in geometry"
		end
		geometry.each_index{|i|
			geometry[i].properties[:charge] = self[i]
		}
		return nil
	end

	def write_file(filename)
		# Write the charges to a text file, format is:
		# number of charges
		# charge type
		# charges, one charge per line
		f = File.open(filename, "w+")
		f.puts size
		f.puts @charge_type
		each{|ch|
			f.puts ch
		}
		f.close
	end

	def write_file_xyzc(filename, geometry)
		# Write the charges to a text file, format is:
		# number of charges
		# charge type
		# x,y,z,charge
		f = File.open(filename, "w+")
		f.puts size
		f.puts @charge_type
		each_index{|i|
			ch = at(i)
			x,y,z = geometry[i].to_a
			f.puts sprintf("%16.8f%16.8f%16.8f%12.6f",x,y,z,ch)
		}
		f.close
	end

	def AtomicCharges.from_file(filename)
		f = File.open(filename, "r")
		begin
			chsize = f.gets
		rescue
			Cuby::error "Error reading atomic charges file '#{filename}'\nnumber of charges not found"
		end
		chsize = chsize.to_i
		if chsize == 0
			Cuby::error "Error reading atomic charges file '#{filename}'\nno information on number of charges in file, or the size is 0"
		end
		begin
			type = f.gets
		rescue
			Cuby::error "Error reading atomic charges file '#{filename}'\ncharge type not found"
		end
		charges = AtomicCharges.new(type.strip.to_sym)
		chsize.times{|i|
			begin
				ch = f.gets
			rescue
				Cuby::error "Error reading atomic charges file '#{filename}'\nless charges than specified in header"
			end
			charges << ch.to_f
		}

		f.close

		return charges
	end

	def AtomicCharges.from_geometry(geometry)
		# read charges from geometry, from Atom.properties[:charge]
		charges = AtomicCharges.new(:unknown)
		geometry.each_with_index{|atom, i|
			unless ch = atom.properties[:charge]
				Cuby::error "Requested using atoimic charges from geometry but charge on atom #{i+1} is not set"
			end
			charges[i] = ch
		}
		return charges
	end
end

class AtomicChargesGradient < Array

	def initialize(n)
		# Upon initialization, zero gradient is built for system of n atoms
		n.times{|i|
			x = []
			n.times{|j| x << Coordinate.new(0,0,0)}
			self[i] = x
		}
	end

	def to_s
		# Nice conversion to string for printing
		s = ""
		size.times{|i|
			s += "atom #{i+1}:\n"
			size.times{|j| 
				s += sprintf("%15.8f", self[i][j].x)
				s += sprintf("%15.8f", self[i][j].y)
				s += sprintf("%15.8f", self[i][j].z)
				s += "\n"
			}
		}
		return s
	end

	def *(arg)
		# Multiplication
		result = AtomicChargesGradient.new(self.size)
		each_index{|i|
			each_index{|j|
				result[i][j] = self[i][j] * arg
			}
		}
		return result
	end
end

