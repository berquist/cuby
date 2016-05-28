require "classes/math/polynomial_from_conditions"

class SoftCoreRepulsion
	# Curve built from three segments:
	# exponential at long distances (above r_exp)
	# constant at short distances (below r_plateau)
	# smooth transition in between

	# the exponential is k * Math::exp(exp * r)

	def initialize(r_plateau, r_exp, k, exp, exp2 = 1.0)
		@r_plateau, @r_exp, @k, @exp, @exp2 = r_plateau, r_exp, k, exp, exp2

		# Transition to exponential
		eex = k * Math::exp(exp * @r_exp**@exp2)
		dexp = @exp * @exp2 * @k * @r_exp**(@exp2-1) * Math::exp(@exp * @r_exp**@exp2)
		d2exp = @exp**2 * @exp2**2 * @k * @r_exp**(2*@exp2-2) * Math::exp(@exp*@r_exp**@exp2) +
			@exp*(@exp2-1) * @exp2 * @k * @r_exp**(@exp2-2) * Math::exp(@exp*@r_exp**@exp2)

		# Value at the end of the plateau
		@e0 = k * Math::exp(exp * @r_plateau**@exp2)

		@transition = PolynomialFromConditions.build([
			PolynomialFromConditions::Condition.new(:value, @r_plateau, @e0),
			PolynomialFromConditions::Condition.new(:first_deriv, @r_plateau, 0.0),
			PolynomialFromConditions::Condition.new(:second_deriv, @r_plateau, 0.0),
			PolynomialFromConditions::Condition.new(:value, @r_exp, eex),
			PolynomialFromConditions::Condition.new(:first_deriv, @r_exp, dexp),
			PolynomialFromConditions::Condition.new(:second_deriv, @r_exp, d2exp)
		])

		@transition_d = @transition.derivative
	end

	def calculate(r)
		if r <= @r_plateau
			return @e0
		elsif r <= @r_exp
			return @transition.calc(r)
		else
			return @k * Math::exp(@exp * r**@exp2)
		end
	end

	def deriv(r)
		if r <= @r_plateau
			return 0.0
		elsif r <= @r_exp
			return @transition_d.calc(r)
		else
			return  @exp * @exp2 * @k * r**(@exp2-1) * Math::exp(@exp * r**@exp2)
		end
	end

	def plot(range = 5.0, step = 0.05)
		x = 0.0
		while x <= range
			printf("%10.5f%10.5f\n", x, calculate(x))
			x += step
		end
	end
end

