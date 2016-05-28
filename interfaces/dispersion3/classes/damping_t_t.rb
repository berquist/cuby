# Tang and Toennies damping function

require "classes/math/factorial.rb"

module DampingTT

	def DampingTT.f6(r,b)
		sum = 0.0
		(0..6).each{|k|
			sum += (b * r)**k / Factorial.n(k)
		}

		return 1.0 - Math::exp(-b * r) * sum
	end

	def DampingTT.f8(r,b)
		sum = 0.0
		(0..8).each{|k|
			sum += (b * r)**k / Factorial.n(k)
		}

		return 1.0 - Math::exp(-b * r) * sum
	end

	def DampingTT.sample_plot(b, rmax = 5.0, step = 0.1)
		r = 0.0
		while r <= rmax do
			printf("%10.3f%10.3e\n", r, DampingTT.f6(r,b))
			r += step
		end
	end
end

#DampingTT.sample_plot(10.0)
