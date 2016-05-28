require "classes/objects/hessian.rb"

class TranslationRotationRemoval

	def initialize(settings)
		#: The settings are read in order to determine what degrees of freedom should
		#: be frozen.
		if settings[:remove_translation] && settings[:remove_rotation]
			@task = :trans_rot
		elsif !settings[:remove_translation] && settings[:remove_rotation]
			@task = :rot
		elsif settings[:remove_translation] && !settings[:remove_rotation]
			@task = :trans
		else
			@task = :none
		end
	end

	def fix_gradient!(geometry, gradient)
		#: Apply the desired removal opration on Gradient object
		case @task
		when :trans_rot
			#: Translation and rotation are removed at the same time
			#: using a projector matrix
			projector = Hessian.transrot_removal_projector(geometry)
			g = projector * gradient.to_vector
			gradient.update_from_vector!(g)
		when :rot
			#: Rotation are removed using a projector matrix
			projector = Hessian.rot_removal_projector(geometry)
			g = projector * gradient.to_vector
			gradient.update_from_vector!(g)
		when :trans
			# Translation removal is simple, projector is not needed
			remove_translation!(gradient)
		end
		return nil
	end

	def remove_translation!(gradient)
		total_grad = Coordinate.new
		gradient.each{|atom|
			total_grad.plus! atom
		}
		gradient.each{|atom|
			atom.minus! total_grad / gradient.size
		}
	end

end
