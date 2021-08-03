using ComponentArrays
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

const G = 6.67408e-11 * 1e-9

struct Body
    mass::Float64                   # [kg]
    radius::Float64                 # [km]
    position::SVector{3, Float64}   # [km]
    speed::SVector{3, Float64}      # [km·s⁻¹]
end

function create_solar_system()
    sun = Body(
        1.989e30,
        695500.,
        SA[0., 0., 0.],
        SA[0., 0., 0.]
    )

    earth = Body(
        5.972e24,
        63710.,
        SA[149_598_023., 0., 0.],
        SA[0., 29.78, 0.]
    )

    moon = Body(
        7.347e22,
        1736.,
        earth.position + SA[0., 406_300., 0],
        earth.speed + SA[-1.022, 0., 0.]
    )
    [sun, earth, moon], [:sun, :earth, :moon]
end


function step!(du, u, p, t)
    masses, names = p

    for body in names
        @view(du[body]).x = u[body].ẋ
        @view(du[body]).ẋ .= 0.

    end
    
    for (i, body₁) in enumerate(names)
        for body₂ in names[i+1:end]
            u_1_2 = u[body₂].x - u[body₁].x
            f_1_2 = G / norm(u_1_2)^2 * normalize(u_1_2)

            @view(du[body₁]).ẋ += f_1_2 * masses[body₂].m
            @view(du[body₂]).ẋ -= f_1_2 * masses[body₁].m
        end
    end
end


function simulate()
    tspan = (0., 100_000_000.)
    dt = 1_000.

    (sun, earth, moon), body_names = create_solar_system()

    u₀ = ComponentArray(
        sun=ComponentArray(
            x=sun.position,
            ẋ=sun.speed
        ),
        earth=ComponentArray(
            x=earth.position,
            ẋ=earth.speed
        ),
        moon=ComponentArray(
            x=moon.position,
            ẋ=moon.speed
        )
    )

    params = ComponentArray(
        sun=ComponentArray(
            m=sun.mass
        ),
        earth=ComponentArray(
            m=earth.mass
        ),
        moon=ComponentArray(
            m=moon.mass
        )
    ), body_names
    ode_prob = ODEProblem(step!, u₀, tspan, params)
    solve(ode_prob, Tsit5(), dt=dt, adaptive=false)
end


@userplot SolarSystemPlot
@recipe function f(bp::SolarSystemPlot)
    x, y, i, color = bp.args
    inds = 1:i
    last_imp_inds = 1000

    alphas = vcat(repeat([0.1], max(0, i-last_imp_inds)), repeat([1.], min(i, last_imp_inds)))
    
    @assert length(alphas) == i

    marker --> :circle
    markersize --> Int[zeros(i-1); 5]
    linecolor --> color
    linewidth --> (alphas .* 2 .|> ceil .|> Int)
    linealpha --> alphas
    # xlims --> [-1, 1] * 150_000_000
    # ylims --> [-1, 1] * 150_000_000
    showaxis --> false
    ticks --> false
    xguide --> ""
    yguide --> ""

    x[inds], y[inds]
end


function draw_animation(ode_sol, body_names)
    x_1 = [p[:earth].x[1] for p in ode_sol.u]
    y_1 = [p[:earth].x[2] for p in ode_sol.u]

    x_2 = [p[:moon].x[1] for p in ode_sol.u]
    y_2 = [p[:moon].x[2] for p in ode_sol.u]
    anim = @animate for i in 100:20:10_000
        solarsystemplot(x_1, y_1, i, :blue)
        solarsystemplot!(x_2, y_2, i, :red)
        xlims!((x_1[i] - 5000_000, x_1[i] + 500_0000))
        ylims!((y_1[i] - 5000_000, y_1[i] + 500_0000))
    end every 5

    gif(anim, "solar.gif", fps=15)
end