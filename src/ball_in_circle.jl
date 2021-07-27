using ComponentArrays
using DifferentialEquations
using LinearAlgebra
using Unitful
using UnitfulRecipes
using Plots

r = 1.0u"m"
g = 9.81u"m/s^2"
t = 30u"s"
dt = 3e-3u"s"

function ball!(du, u, p, t)
    # @show u, t
    du.x .= u.ẋ
    du.ẋ[1] = 0.0u"m/s^2"
    du.ẋ[2] = -g
end

function bounce_condition(u, t, integrator) 
    (r^2 - sum(u.x.^2)) / 1u"m^2"
end

function bounce!(integrator)
    x = integrator.u.x
    ẋ = integrator.u.ẋ
    integrator.u.ẋ = ẋ - 2 * (ẋ' * x / norm(x)^2 .* x)
end

function solve_ode(x₀, y₀, ẋ₀, ẏ₀)
    u₀ = ComponentArray(x=[x₀, y₀]u"m", ẋ=[ẋ₀, ẏ₀]u"m/s")
    tspan = (0.0u"s", t)
    
    bounce_cb = ContinuousCallback(bounce_condition, bounce!)
    
    prob = ODEProblem(ball!, u₀, tspan, nothing)
    solve(prob, Tsit5(), callback=bounce_cb, dt=dt, adaptive=false)
end

sol_1  = solve_ode(0.001, 0., 0., 0.)
sol_2  = solve_ode(0.00105, 0., 0., 0.)

@userplot BouncingPlot
@recipe function f(bp::BouncingPlot)
    x, y, i, color = bp.args
    n = 1u"s" / 4dt |> floor |> Int

    if i > n
        inds = i - n + 1:i |> collect
    else
        inds = vcat(repeat([1], n - i), 1:i)
    end

    marker --> :circle
    markercolor --> color
    markerstrokecolor --> color
    markersize --> 3
    seriesalpha --> range(0.2, 1, length=n)
    xlims --> [-1, 1]u"m"
    ylims --> [-1, 1]u"m"
    showaxis --> false
    ticks --> false
    xguide --> ""
    yguide --> ""

    x[inds], y[inds]
end

function circle_shape(h, k, r)
    θ = range(0, 2π, length=500)
    h .+ r * sin.(θ), k .+ r * cos.(θ)
end

n = t / dt |> floor |> Int
x_1 = [uu.x[1] for uu in sol_1.u]
y_1 = [uu.x[2] for uu in sol_1.u]

x_2 = [uu.x[1] for uu in sol_2.u]
y_2 = [uu.x[2] for uu in sol_2.u]

anim = @animate for i in 1:n
    bouncingplot(x_1, y_1, i, :blue)
    bouncingplot!(x_2, y_2, i, :red)
    plot!(circle_shape(0u"m", 0u"m", 1u"m"), seriestype=[:shape], linecolor=:black, legend=false, aspect_ratio=1, fillalpha=0.1, c=:blue)
end every 3

gif(anim, "unitful_anim.gif", fps=30)