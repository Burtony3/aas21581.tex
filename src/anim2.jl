### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 832a3aa9-3608-439f-a71f-c612852d0a93
using Javis

# ╔═╡ 1ed3136f-e2a6-4609-b4ab-91f7e9e834fa
function object(r, p = O, color = "black")
    sethue(color)
    circle(p, r, :fill)
    return p
end

# ╔═╡ a7f6a3b7-0ddc-424e-b1d2-c4a8a23aa200
function connector!(connection, p1, p2, color)
    sethue(color)
    push!(connection, [p1, p2])
    map(x -> line(x[1], x[2], :stroke), connection)
end

# ╔═╡ 6ece23ff-beb1-4664-91d6-cb018f4ee727
function ground(args...)
    background(40/255, 44/255, 52/255) # canvas background
    sethue(204/255, 204/255, 204/255) # pen color
end

# ╔═╡ 799b768e-e4e6-47a8-9afa-b8091721613b
function circ(p = O, color = "black", action = :fill, radius = 25, edge = "solid")
    sethue(color)
    setdash(edge)
    circle(p, radius, action)
end

# ╔═╡ a27dd57d-c34f-47a2-874b-70a5db688a96
begin

    # to store the connectors
    connection = []

    frames = 1000

    # setup the video
    myvideo = Video(500, 500)
    Background(1:frames, ground)
	
	# colors
	cred = (191/255, 97/255, 106/255)
	cwhite = (204/255, 204/255, 204/255)
	cblue = (129/255, 161/255, 193/255)
	cpurple = "#b48ead"

    # draw the orbits
    earth_orbit = Object((args...) -> circ(O, cwhite, :stroke, 200))
    venus_orbit = Object((args...) -> circ(O, cwhite, :stroke, 144))

    # add the objects
    earth = Object(1:frames, (args...) -> object(5, O, cblue), Point(200, 0))
    venus = Object(1:frames, (args...) -> object(4, O, cred), Point(144, 0))

    # move the planets
    # We need the planets to revolve according to their time periods.
    # Earth completes its one revolution in 365 days and Venus does that in 224.7 days.
    # Hence, we need to multiply (224.7/365) so that the time period matches properly i.e.,
    # when earth completes its full revolution, Venus has done (224.7/365) th of its revolution.
    act!(earth, Action(anim_rotate_around(12 * 2π * (224.7 / 365), O)))
    act!(venus, Action(anim_rotate_around(12 * 2π, O)))

    # draw the connectors
    Object(1:frames, (args...) -> connector!(connection, pos(earth), pos(venus), cpurple))

    # render
    render(myvideo, pathname = "cosmic_dance2.gif")
end

# ╔═╡ Cell order:
# ╠═832a3aa9-3608-439f-a71f-c612852d0a93
# ╠═1ed3136f-e2a6-4609-b4ab-91f7e9e834fa
# ╠═a7f6a3b7-0ddc-424e-b1d2-c4a8a23aa200
# ╠═6ece23ff-beb1-4664-91d6-cb018f4ee727
# ╠═799b768e-e4e6-47a8-9afa-b8091721613b
# ╠═a27dd57d-c34f-47a2-874b-70a5db688a96
