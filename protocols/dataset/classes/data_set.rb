require "protocols/dataset/classes/plot.rb"

#=======================================================================
# Dataset definition classes, used in the yaml representation of the set
#=======================================================================

module ProtocolDataset
	class DataSetDescription
		attr_reader :name
		attr_reader :dataset_setup # hash, to be copied into the settings
		attr_reader :global_setup # hash, to be copied into calculation settings
		attr_reader :groups
		attr_reader :references # hash DOI -> reference
		attr_reader :text

		def initialize
			# Sample only, real data are read from yaml
			@name = "Data set example"
			@dataset_setup = {}
			@global_setup = {}
			@groups = ["group"] # is nil where there are no groups
		end
	end

	class DataSetItem
		attr_reader :name
		attr_reader :shortname
		attr_reader :geometry # file or yaml
		attr_accessor :reference_value
		attr_reader :setup # hash, to be copied into calculation settings
		attr_reader :group
		attr_reader :tags # string
		attr_accessor :add_to_result

		def initialize
			# Sample only, real data are read from yaml
			@name = "System example"
			@shortname = "dataset_example"
			@geometry = "geofile.xyz"
			@reference_value = 0.0
			@setup = {}
			@group = "group"
			@tags = "example, item"
			@add_to_result = nil
		end
	end

	class DataSetPlot
		# A plot (e.g. dissociation curve) in dataset
		attr_reader :name # String, printed in the plot
		attr_reader :filename # String, filename-friendly
		attr_reader :first # Index of first item in the series
		attr_reader :last # Index of last item in the series
		attr_reader :x_values # Array of the x values

		def save(dataset, jobs, format, one_file, extra_cols, methodname)
			x = []
			y_ref = []
			y = []
			@x_values.each_index{|i|
				item_i = @first + i
				x << @x_values[i]
				y_ref << dataset.original_items[item_i].reference_value
				y << jobs[dataset.items.index(dataset.original_items[item_i])].results.energy
			}
			case format
			when :text
				f = File.open(@filename+".txt", "w+")
				@x_values.each_index{|i|
					f.printf("%12.6f%12.6f%12.6f\n", x[i], y_ref[i], y[i])
				}
				f.close
			when :gnuplot
				y_cols = [y_ref, y]
				y_names = ["reference", methodname]
				extra_cols.each{|series|
					y_names << series[0]
					y_cols << series[1]
				}
				plot = Plot.new(x, y_cols, y_names)
				plot.title = @name
				plot.gnuplot_png(@filename)
			when :yaml
				one_file.puts ({@name => y}).to_yaml.gsub(/^---$/,"")
			else
				raise "Unknown value of plot format"
			end
		end

		def is_complete?(dataset)
			@x_values.each_index{|i|
				item_i = @first + i
				return false unless dataset.items.index(dataset.original_items[item_i])
			}
			return true
		end
	end

	class DataSet
		attr_reader :description
		attr_accessor :items
		attr_accessor :alternative_reference # Hash of Arrays
		attr_accessor :plots
		attr_reader :original_items # Preserves the whole set before a selection is made

		def initialize
			# Sample only, real data are read from yaml
			@description = DataSetDescription.new
			@items = [DataSetItem.new]
			@alternative_reference = nil
		end

		def save_orig_items
			@original_items = []
			@items.each_index{|i|
				@original_items[i] = @items[i]
			}
			return nil
		end

		def save_plots(jobs, settings)
			return unless @plots
			format = settings[:dataset_save_plots]
			format = :gnuplot if settings[:dataset_save_plots] == :gnuplot_tiled

			# Get plots to be built
			active_plots = []
			@plots.each{|plot|
				# Check if the plot is complete in this selection
				active_plots << plot if plot.is_complete?(self)
			}
			if active_plots.empty?
				Cuby::warning("The selection applied does not yield complete curves for plotting")
				return
			end

			# Add extra cols from input
			extra_cols = {}
			active_plots.each{|plot| extra_cols[plot.name] = []}

			settings[:dataset_plot_extra].each_pair{|name, filename|
				f = File.open(filename)
				plotdata = YAML.load(f)
				f.close

				active_plots.each{|plot|
					extra_cols[plot.name] << [name, plotdata[plot.name]]
				}
			}

			# Writting all into one file - prepare the file
			if settings[:dataset_save_plots] == :yaml
				one_file = File.open(settings[:dataset_plots_fn] + ".yaml", "w+")
			else
				one_file = nil
			end

		
			# Save plots
			active_plots.each_index{|i|
				plot = active_plots[i]
				plot.save(self, jobs, format, one_file, extra_cols[plot.name], settings[:dataset_plot_methodname])
			}

			# Writting all into one file - close the file
			one_file.close if one_file

			# Tile images using imagemagick
			if settings[:dataset_save_plots] == :gnuplot_tiled
				filelist = ""
				active_plots.each{|plot| filelist += plot.filename + ".png " }

				begin
					# Tile images
					`montage #{filelist} -geometry +2+2 -tile #{settings[:dataset_plot_columns]}x #{settings[:dataset_plots_fn]}.png`
					# delete original files
					`rm -f #{filelist}`
				rescue
					Cuby::error "Calling imagemagick for composing the plots failed"
				end
			end
		end
	end
end
