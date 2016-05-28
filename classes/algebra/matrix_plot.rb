
module MatrixPlot

	#=======================================================================
	# Plot as image using gnuplot
	#=======================================================================

	def plot_gnuplot(filename, run_gnuplot = false, viewer = nil)
		# Write data
		plot_gnuplot_data(filename)
		# Write gnuplot script
		plot_gnuplot_script(filename)
		# Run gnuplot
		if run_gnuplot
			system "gnuplot < #{filename}.gnuplot 2> /dev/null > /dev/null"
		end

		if viewer
			system "#{viewer} #{filename}.png 2> /dev/null > /dev/null"
		end
	end

	def plot_gnuplot_data(filename)
		f = File.open(filename + ".dat","w+")
		m.times{|i|
			f.puts(self.row_as_array(i).join(" "))
		}
		f.close
		return nil
	end

	def plot_gnuplot_script(filename)
		f = File.open(filename + ".gnuplot","w+")
		f.puts "reset"
		f.puts "set term png"
		f.puts "set size ratio #{self.m.to_f / self.n}"
		f.puts "set xrange [-1:#{self.n}]"
		f.puts "set yrange [-1:#{self.m}] reverse"
		# Minimum, maximum
		min = self.min
		max = self.max
		if min >= 0 && max >= 0
			# Start at zero
			min = 0
		else
			# Symmetric scale
			min = -max
		end

		f.puts "set cbrange [#{min}:#{max}]"
		# RGB palette
		f.puts "set palette defined ( 0 0.0 0.0 1.0, 0.25 0.0 1.0 1.0, 0.5 0.0 1.0 0.0, 0.75 1.0 1.0 0.0, 1.0 1.0 0.0 0.0)"
		f.puts "unset key"
		f.puts "set out '#{filename}.png'"
		f.puts "plot '#{filename}.dat' matrix with image"
		f.close
		return nil
	end

end
