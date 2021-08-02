using ComponentArrays
using DifferentialEquations
using LinearAlgebra
using Plots
using StaticArrays

const G = 6.67408e-11 * 1e-9

mutable struct Body
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
    [sun, earth]
end


function step!(du, u, p, t)
    du.sun.x = u.sun.ẋ
    du.earth.x = u.earth.ẋ

    u_sun_earth = u.earth.x - u.sun.x
    f_sun_earth = G / norm(u_sun_earth)^2 * normalize(u_sun_earth)

    du.sun.ẋ = f_sun_earth * p.earth.m
    du.earth.ẋ = - f_sun_earth * p.sun.m

    @show u, t
end


function simulate()
    tspan = (0., 100_000.)
    dt = 1000.

    sun, earth = create_solar_system()

    u₀ = ComponentArray(
        sun=ComponentArray(
            x=sun.position,
            ẋ=sun.speed
        ),
        earth=ComponentArray(
            x=earth.position,
            ẋ=earth.speed
        )
    )

    params = ComponentArray(
        sun=ComponentArray(
            m=sun.mass
        ),
        earth=ComponentArray(
            m=earth.mass
        )
    )
    ode_prob = ODEProblem(step!, u₀, tspan, params)
    solve(ode_prob, Tsit5(), dt=dt, adaptive=false)
end