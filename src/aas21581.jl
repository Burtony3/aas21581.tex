module aas21581
using DifferentialEquations, LinearAlgebra, Plots

# ============== #
#   DATA TYPES   #
# ============== #
struct CisLunarTransfer
    Type
    Subtype
    v∞
    θ₀
    λ₀
    dist
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
const tofmax = [40.0, 67.5, 95.0, 120.0, 165.0, 210.0]
const typeDict = Dict(:A => 1, :B => 2, :C => 3, :D => 4, :E => 5, :F => 6)
const subtypeDict = Dict(:ii => (-1, -1), :io => (-1, 1), :oi => (1, -1), :oo => (1, 1))

# =========== #
#   EXPORTS   #
# =========== #
export Trajectory, plot, Optimize, show

# ====================== #
#   INTERNAL FUNCTIONS   #
# ====================== #

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

function cb_condition(u, t, integrator)
    dₘ - norm(u[1:3])
end

function cb_save!(integrator, array)
    push!(array, [integrator.t, integrator.u...])
end

function SetupODE(v∞, θ₀, λ₀, type, subtype, fixed_tof=nothing)
    # INITIALIZING LUNAR STATES UNIT VECTORS
    r̂ₘ = [cosd(θ₀), sind(θ₀), 0]
    v̂ₘ = [-sind(θ₀), cosd(θ₀), 0]

    # CONVERTING INPUTS TO INITIAL CONDITIONS
    outgoing, incoming = subtypeDict[subtype]
    tofidx = typeDict[type]
    r⃗₀ = dₘ*r̂ₘ
    v⃗₀ = v∞*cosd(λ₀)*sind(outgoing*90)*r̂ₘ + (v∞*sind(λ₀) + vₘ)*v̂ₘ
    if isnothing(fixed_tof)
        tof = tofmax[tofidx]*86400.0
    else
        tof = fixed_tof
    end
    tspan = (0.0, tof)

    # CREATING THE PROBLEM
    prob = ODEProblem(CR3BP!, vcat(r⃗₀, v⃗₀), tspan)
end

function GenerateJacobian(traj)

end

# ======================= #
#   EXPORTED FUNCTIONS    #
# ======================= #

function Trajectory(v∞, θ₀, λ₀, type, subtype)

    # SETTTING UP ODE PROBLEM & SOLVING
    x = []
    cb = ContinuousCallback(cb_condition, (integrator) -> cb_save!(integrator, x), nothing, save_positions=(false, false))
    sol = solve(SetupODE(v∞, θ₀, λ₀, type, subtype), reltol=1e-8, abstol=1e-8, callback=cb)

    # SAVING LUNAR CROSSINGS
    x = hcat(x...)
    dist = []
    outgoing, incoming = subtypeDict[subtype]
    tofidx = typeDict[type]
    if ~isempty(x)
        # LOOPING CALLBACK RETURNS
        for i = 1:size(x)[2]
            r⃗ₘ_end = dₘ*[cosd(θ₀ + ωₘ*x[1, i]), sind(θ₀ + ωₘ*x[1, i]), 0]
            d = r⃗ₘ_end - [x[2, i], x[3, i], x[4, i]]
            r⃗ = [x[2, i], x[3, i], x[4, i]]
            v⃗ = [x[5, i], x[6, i], x[7, i]]

            # CHECKING IF WITHIN TOF AND DIRECTION BOUNDS
            if x[1, i] > (tofmax[tofidx-1]-10)*86400# && incoming*(acosd(dot(r⃗, v⃗)/(norm(r⃗)*norm(v⃗)))-90) < 0
                push!(dist, [x[1, i], norm(d)])
            end # if x[1, i]...
        end # for i
    end # if ~isempty
    # APPENDING END POINT DISTANCE
    r⃗ₘ_end = dₘ*[cosd(θ₀ + ωₘ*sol.t[end]), sind(θ₀ + ωₘ*sol.t[end]), 0]
    d = r⃗ₘ_end - sol.u[end][1:3]
    push!(dist, [sol.t[end], norm(d)])

    # OUTPUTTING
    traj = CisLunarTransfer(type, subtype, v∞, θ₀, λ₀, dist, transpose(hcat(sol.u...)), sol.t)
end

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

# ============================ #
#   MODIFYING OTHER PACKAGES   #
# ============================ #

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

function Base.show(io::IO, m::CisLunarTransfer) 
    println("Family: ", String(m.Type), String(m.Subtype))
    println("Launch v∞: ", m.v∞, " km/s || Latitude: ", m.λ₀, " Degrees")
    println("Sun-Earth-Moon Angle: ", m.θ₀, " Degrees")
    min = Inf
    for i = 1:length(m.dist)
        if m.dist[i][2] < min; min = m.dist[i][2]; end
    end
    print("Closest Approach: ", min, " km")
end


end