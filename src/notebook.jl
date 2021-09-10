### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 1c9a7102-b399-4431-9799-474899f5d35d
begin
	include("./aas21581.jl")
	using LinearAlgebra, DifferentialEquations, Plots, PlutoUI
	import .aas21581
	
end

# ╔═╡ f552e7fb-66da-4534-96f8-46970f317dd9
begin
	using Javis
	using Random

	function ground(args...)
		background("white")
		sethue("black")
	end

	function draw_line(p1 = O, p2 = O, color = "black", action = :stroke, edge = "solid")
		sethue(color)
		setdash(edge)
		line(p1, p2, action)
	end

	function circ(p = O, color = "black", action = :fill, radius = 25, edge = "solid")
		sethue(color)
		setdash(edge)
		circle(p, radius, action)
	end

	function info_box(video, object, frame)
		fontsize(12)
		box(140, -210, 170, 40, :stroke)
		Text("10-20 EEG Array Readings", 140, -220, valign = :middle, halign = :center)
		Text("t = $(frame)s", 140, -200, valign = :middle, halign = :center)
	end

	function electrode(
		p = O,
		fill_color = "white",
		outline_color = "black",
		action = :fill,
		radius = 25,
		circ_text = "",
	)
		sethue(fill_color)
		circle(p, radius, :fill)
		sethue(outline_color)
		circle(p, radius, :stroke)
		Text(circ_text, p, valign = :middle, halign = :center)
	end

	electrodes_list = [
		(name = "Cz", position = O),
		(name = "C3", position = Point(-70, 0)),
		(name = "C4", position = Point(70, 0)),
		(name = "T3", position = Point(-140, 0)),
		(name = "T4", position = Point(140, 0)),
		(name = "Pz", position = Point(0, 70)),
		(name = "P3", position = Point(-50, 70)),
		(name = "P4", position = Point(50, 70)),
		(name = "Fz", position = Point(0, -70)),
		(name = "F3", position = Point(-50, -70)),
		(name = "F4", position = Point(50, -70)),
		(name = "F8", position = Point(115, -80)),
		(name = "F7", position = Point(-115, -80)),
		(name = "T6", position = Point(115, 80)),
		(name = "T5", position = Point(-115, 80)),
		(name = "Fp2", position = Point(40, -135)),
		(name = "Fp1", position = Point(-40, -135)),
		(name = "A1", position = Point(-190, -10)),
		(name = "A2", position = Point(190, -10)),
		(name = "O1", position = Point(-40, 135)),
		(name = "O2", position = Point(40, 135)),
	]

	radius = 15
	indicators = ["white", "gold1", "darkolivegreen1", "tomato"]
	demo = Video(500, 500)

	anim_background = Background(1:10, ground)
	head = Object((args...) -> circ(O, "black", :stroke, 170))
	inside_circle = Object((args...) -> circ(O, "black", :stroke, 140, "longdashed"))
	vert_line = Object(
		(args...) ->
			draw_line(Point(0, -170), Point(0, 170), "black", :stroke, "longdashed"),
	)
	horiz_line = Object(
		(args...) ->
			draw_line(Point(-170, 0), Point(170, 0), "black", :stroke, "longdashed"),
	)

	for num in 1:length(electrodes_list)
		Object(
			(args...) ->
				electrode.(
					electrodes_list[num].position,
					rand(indicators, length(electrodes_list)),
					"black",
					:fill,
					radius,
					electrodes_list[num].name,
				),
		)
	end
	info = Object(info_box)

	render(demo, pathname = "test.gif", framerate = 1)
end

# ╔═╡ b08cb836-7c51-4b92-9269-0e738ef1db24
T = aas21581.Trajectory(0.9, 0, -12.9, :C, :oi)

# ╔═╡ 22e3e298-c45b-4cbb-acad-33713b4a8aa2
begin
	aas21581.plot(T, 0)
	# legend("off")
	# xaxis!("X-Position (km) - Inertial Earth Centered")
	# yaxis!("Y-Position (km)")
	# png("C:/Users/Burto/repos/aas21581.tex/etc/oo.png")
end

# ╔═╡ 5aedf9df-5e82-4627-b5c8-26a9a68cc385
T2 = aas21581.Optimize(T)

# ╔═╡ 0ac1c328-8e87-4dcd-8acf-cb1c2a10b47a
aas21581.plot(T2);

# ╔═╡ Cell order:
# ╠═1c9a7102-b399-4431-9799-474899f5d35d
# ╠═b08cb836-7c51-4b92-9269-0e738ef1db24
# ╠═22e3e298-c45b-4cbb-acad-33713b4a8aa2
# ╠═5aedf9df-5e82-4627-b5c8-26a9a68cc385
# ╠═0ac1c328-8e87-4dcd-8acf-cb1c2a10b47a
# ╠═f552e7fb-66da-4534-96f8-46970f317dd9
