
module DefaultInterfaceSelection
	@@methods_printed = []

	def DefaultInterfaceSelection.from_settings(settings)
		if settings[:default_interfaces][settings[:method]]
			interface = settings[:default_interfaces][settings[:method]]
			unless @@methods_printed.include?(settings[:method])
				Cuby::log.puts_v(:normal, "Default interface #{interface} used for method #{settings[:method]}")
				@@methods_printed << settings[:method]
			end
			return interface
		else
			Cuby::error("Interface can not be determined automatically,\nno default interface set for method #{settings[:method]}")
		end
	end
end
