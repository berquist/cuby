class Plot
	attr_reader :x
	attr_reader :y_columns
	attr_reader :y_names
	attr_accessor :x_label
	attr_accessor :y_label
	attr_accessor :title
	attr_accessor :x_range # optional, array of two numbers - from and to
	attr_accessor :y_range # optional, array of two numbers - from and to
	attr_accessor :zero_line 

	#=======================================================================
	# Class setup
	#=======================================================================
	
	LINE_STYLES_DEFAULT = [
		"set style line 1 lc rgb 'black'",
		"set style line 2 lc rgb 'red'",
		"set style line 3 lc rgb 'blue'",
		"set style line 4 lc rgb 'green'"
	]

	@@line_styles = LINE_STYLES_DEFAULT

	def Plot.set_line_styles(styles)
		@@line_styles = styles
	end

	#=======================================================================
	# methods
	#=======================================================================
	
	def initialize(x, y_columns, y_names)
		@x = x
		@y_columns = y_columns
		@y_names = y_names
	end

	def gnuplot_script(file, data_file_name)
		# Line width
		lw = 2.0

		#file.puts "set term png"
		file.puts "set terminal pngcairo dashed"
		file.puts "set termoption dash"
		file.puts "set title '#{@title}'" if @title
		file.puts "set xlabel '#{@x_label}'" if @x_label
		file.puts "set ylabel '#{@y_label}'" if @y_label
		file.puts "set xrange [#{@x_range[0]}:#{@x_range[1]}]" if @x_range
		file.puts "set yrange [#{@y_range[0]}:#{@y_range[1]}]" if @y_range

		# Colors and styles
		@@line_styles.each{|style| file.puts style}

		# zero line
		if @zero_line
			raise "x_range must be set when zero_line is used" unless @x_range
			file.puts "set style line 100 lt 1 lc rgb 'gray' lw #{lw}"
			file.puts "set arrow 100 from #{@x_range[0]},0.0 to #{@x_range[1]},0.0 nohead ls 100 back"
		end

		file.puts "plot \\"
		lines = []
		y_names = @y_names # Some escaping could be added here
		@y_columns.each_index{|i|
			lines << "'#{data_file_name}' using 1:#{i+2} title \"#{y_names[i]}\" w linespoints ls #{i+1}"
		}
		file.puts lines.join(",\\\n")

	end

	def gnuplot_data_file(file)
		@x.each_index{|i|
			file.printf("%12.6f", @x[i])
			y_columns.each_index{|ci|
				file.printf("%12.6f", @y_columns[ci][i])
			}
			file.puts unless i == @x.size - 1
		}
	end

	def gnuplot_png(filename)
		tmp_fn = "temp_#{Time.now.to_f}"
		# Write data
		file = File.open(tmp_fn + ".dat", "w+")
		gnuplot_data_file(file)
		file.close
		# Write gnuplot script
		file = File.open(tmp_fn + ".gnuplot", "w+")
		file.puts "set out '#{filename}.png'"
		gnuplot_script(file, tmp_fn + ".dat")
		file.close
		# Run gnuplot
		#system "gnuplot < #{tmp_fn}.gnuplot 2> /dev/null > /dev/null"
		system "gnuplot < #{tmp_fn}.gnuplot "
		# Delete temp files
		File.delete(tmp_fn + ".dat")
		File.delete(tmp_fn + ".gnuplot")
	end
end
