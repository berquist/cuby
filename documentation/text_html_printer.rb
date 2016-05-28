
class TextHtmlPrinter
	attr_reader :html

	def initialize(file = $stdout, format = :text)
		raise "Unknown format, use :text or :html" unless [:text, :html].include? format
		@f = file
		@html = format == :html
	end

	def header
		if @html
			#@f.puts "<html>"
			#@f.puts "<head>"
			#@f.puts "</head>"
			#@f.puts "<body>"
		end
	end

	def footer
		if @html
			#@f.puts "</body>"
			#@f.puts "</html>"
		else
			@f.puts
		end
	end

	def heading(text)
		if @html
			@f.puts "<h1>#{text}</h1>"
		else
			@f.puts "================================================================================"
			@f.puts " #{text}"
			@f.puts "================================================================================"
			@f.puts
		end
	end

	def caption(text)
		if @html
			@f.puts "<h2>#{text}</h2>"
		else
			@f.puts "--------------------------------------------------------------------------------"
			@f.puts " #{text}"
			@f.puts "--------------------------------------------------------------------------------"
			@f.puts
		end
	end

	def paragraph(text)
		if @html
			@f.puts "<p>#{text}</p>"
		else
			@f.puts "#{text}"
			@f.puts
		end
	end

	def list(title, array)
		if @html
			@f.puts "<h4>#{title}</h4>" if title
			@f.puts "<ul>"
			array.each{|item| @f.puts "<li>#{item}</li>"}
			@f.puts "</ul>"
		else
			@f.puts "#{title}" if title
			array.each{|item| @f.puts "   * #{item}"}
			@f.puts
		end
	end
end
