### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 7d9fc98e-ea4c-11eb-3478-d156c825a95b
using DifferentialEquations, LinearAlgebra, Plots

# ╔═╡ f659428a-cde8-412a-97b6-e7e20e89b701
begin
	# ============== #
	#   DATA TYPES   #
	# ============== #
	struct CisLunarTransfer
		Type
		Subtype
		v∞
		θ₀
		λ₀
		Crossings
		States
		Times
	end

	# ============= #
	#   CONSTANTS   #
	# ============= #
	const μₑ = 3.986592936294783e5
	const μₛ = 1.327124400419400e11
	const dₛ = 150.0e6
	const ωₛ = sqrt((μₑ + μₛ)/dₛ^3)*180/π
	const μₘ = 4.843941639988467e3
	const dₘ = 3.844e5
	const vₘ = sqrt((μₑ + μₘ)/dₘ)
	const ωₘ = sqrt((μₑ + μₘ)/dₘ^3)*180/π
	const toflib = [40.0, 67.5, 95.0, 120.0, 165.0, 210.0]
	const type2idx = Dict(:A => 1, :B => 2, :C => 3, :D => 4, :E => 5, :F => 6)
	const subtype2dir = Dict(:ii => (-1, -1), :io => (-1, 1), :oi => (1, -1), :oo => (1, 1))
end

# ╔═╡ 6c31c956-7ccc-4a1a-a421-c34b2a420b9a
function CR3BP!(du, u, p, t)
    # DECODING INPUTS
    r⃗ = u[1:3]
    v⃗ = u[4:6]

    # CALCULATING INTERMEDIATES
    r⃗ₛ = dₛ*[-cosd(ωₛ*t), -sind(ωₛ*t), 0]

    # OUTPUTTING
    du[1:3] = v⃗
    du[4:6] = -μₑ*(r⃗/norm(r⃗)^3) - μₛ*( ((r⃗ - r⃗ₛ)/norm(r⃗ - r⃗ₛ)^3) + (r⃗ₛ/norm(r⃗ₛ)^3) )
end

# ╔═╡ b3135ffd-68d5-474a-948c-d2a7c6368945
begin
	function condition(u, t, integrator)
		dₘ - norm(u[1:3])
	end

	function action!(integrator, tcross, ucross)
		# push!(array, [integrator.t, integrator.u])
		# display(array)
		# vcat
		push!(tcross, integrator.t)
		push!(ucross, integrator.u...)
		# display(integrator.u)
		# display(ucross)
	end
end

# ╔═╡ d6192ac7-818b-40c5-8ea4-5f92f6fd2140
function θccw(x, y, ref=[0, 0, 1])
	θ = atand(norm(cross(x, y)), dot(x, y))
	
	if dot(cross(x, y), ref) < 0; θ = 360-θ; end
	
	return θ
	
end

# ╔═╡ 98b59c43-2756-4d98-9103-a207f9d50e6f
function distance2moon(u, t, θ₀)
	# FINDING DISTANCE
	r⃗ₘ = dₘ*[cosd(θ₀+ωₘ*t), sind(θ₀+ωₘ*t), 0]
	dist = norm(r⃗ₘ - u[1:3])
	
	# FINDING DIRECTION
	ψ = θccw(u[1:3], u[4:6])
	dir = ψ < 270 && ψ > 90 ? -1 : 1
	
	return dist, dir
	
end

# ╔═╡ 89c46619-e8ae-4012-884d-89818add139d
function Plots.plot(traj::CisLunarTransfer, subtidx=0)
    θ = LinRange(0, 360, 100)
	plot(3.844e5*cosd.(θ), 3.844e5*sind.(θ), label="Moon", alpha=0.25, color=:black)
	plot!(traj.States[1:end-subtidx, 1], traj.States[1:end-subtidx, 2], 
		aspect_ratio=:equal, label="Trajectory", 
        color=:green, dpi=300)
    ω₁ = sqrt((μₑ + μₘ)/dₘ^3)*180/π
    r⃗ₘ_end = dₘ*[cosd(traj.θ₀ + ω₁*traj.Times[end]), sind(traj.θ₀ + ω₁*traj.Times[end]), 0]
    # scatter!([traj.States[1, 1], r⃗ₘ_end[1]], [traj.States[1, 2], r⃗ₘ_end[2]], label="", color=:black)
    scatter!([traj.States[1, 1]], [traj.States[1, 2]], label="", color=:dodgerblue)
    scatter!([traj.States[end-subtidx, 1]], [traj.States[end-subtidx, 2]], label="", color=:crimson, legend=false, axis=nothing, showaxis=false, fmt=:png)
end

# ╔═╡ bee40bea-5a4c-4160-bb73-0a166f0fa64d
function Trajectory(v∞, θ₀, λ₀, type, subtype)
	println("\n\n==========================")
    # INITIALIZING LUNAR STATES UNIT VECTORS
    r̂ₘ = [cosd(θ₀), sind(θ₀), 0]
    v̂ₘ = [-sind(θ₀), cosd(θ₀), 0]

    # CONVERTING INPUTS TO INITIAL CONDITIONS
    outgoing, incoming = subtype2dir[subtype]
    tofidx = type2idx[type]
    r⃗₀ = dₘ*r̂ₘ
    v⃗₀ = v∞*cosd(λ₀)*sind(outgoing*90)*r̂ₘ + (v∞*sind(λ₀) + vₘ)*v̂ₘ
    tspan = (0.0, toflib[tofidx]*86400.0)

    # CREATING THE PROBLEM
    prob = ODEProblem(CR3BP!, vcat(r⃗₀, v⃗₀), tspan)

    # SETTTING UP ODE PROBLEM & SOLVING
    tcross = []
	ucross = []
    cb = ContinuousCallback(
		condition, 
		(integrator) -> action!(integrator, tcross, ucross), nothing, 
		save_positions=(false, false)
	)
    sol = solve(
		prob, 
		reltol=1e-8, abstol=1e-8, 
		callback=cb
	)
	
	# HANDLING CROSS CONDITIONS
	ucross = reshape(ucross, 6, :)
	distmin = 5e4
	crossing = []
	for i = 1:length(tcross)
		dist, dir = distance2moon(ucross[:, i], tcross[i], θ₀)
		if dist < distmin && dir == incoming && tcross[i]/86400 > toflib[tofidx-1]
			distmin = dist
			crossing = [tcross[i]/86400, dist]
		end
	end
	dist, dir = distance2moon(sol.u[end], sol.t[end], θ₀)
	if dist < distmin && dir == incoming
		crossing = [sol.t/86400, dist]
	end

    # OUTPUTTING
	print("\n")
    traj = CisLunarTransfer(type, subtype, v∞, θ₀, λ₀, crossing, transpose(hcat(sol.u...)), sol.t)
end

# ╔═╡ 4fdc6a36-d972-431e-b6ed-206038cb795c
function Optimize(traj::CisLunarTransfer)
    tof = -1
    for i = 1:length(traj.dist)
        if traj.dist[i][2] < 5e4
            tof = traj.dist[i][1]
        end
    end

    if tof == -1; return nothing; end

    sol = solve(SetupODE(traj.v∞, traj.θ₀, traj.λ₀, traj.Type, traj.Subtype, tof), reltol=1e-8, abstol=1e-8)

    # FINDING FIRST DISTANCE
    ω₁ = sqrt((μₑ + μₘ)/dₘ^3)*180/π
    r⃗ₘ_end = (t) -> dₘ*[cosd(traj.θ₀ + ω₁*t), sind(traj.θ₀ + ω₁*t), 0]
    d = [sol.u[end][1:3]-r⃗ₘ_end(sol.t[end])]
    dist = []
    dist = push!(dist, [tof, norm(d)])

    # while norm(dist) > 500
        

    # end

    traj = CisLunarTransfer(traj.Type, traj.Subtype, traj.v∞, traj.θ₀, traj.λ₀, dist, transpose(hcat(sol.u...)), sol.t)
end

# ╔═╡ e13f4913-3da5-4f40-b092-00dd937830d2
T = Trajectory(0.9, 0, -12.4, :C, :oi)

# ╔═╡ 64e1841d-442b-4ebf-b4ab-ba1ca84a6e87
plot(T)

# ╔═╡ Cell order:
# ╠═7d9fc98e-ea4c-11eb-3478-d156c825a95b
# ╠═f659428a-cde8-412a-97b6-e7e20e89b701
# ╟─6c31c956-7ccc-4a1a-a421-c34b2a420b9a
# ╟─b3135ffd-68d5-474a-948c-d2a7c6368945
# ╟─d6192ac7-818b-40c5-8ea4-5f92f6fd2140
# ╟─98b59c43-2756-4d98-9103-a207f9d50e6f
# ╠═89c46619-e8ae-4012-884d-89818add139d
# ╠═bee40bea-5a4c-4160-bb73-0a166f0fa64d
# ╟─4fdc6a36-d972-431e-b6ed-206038cb795c
# ╠═e13f4913-3da5-4f40-b092-00dd937830d2
# ╠═64e1841d-442b-4ebf-b4ab-ba1ca84a6e87
