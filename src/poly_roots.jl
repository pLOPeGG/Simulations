import Base: +, isless, length

struct Polynom{T <: Number}
    coeffs::AbstractVector{T}
end

function normalize(p::Polynom{T}) where T
    @assert length(p) > 0 && p.coeffs[end] != zero(T)
    Polynom(p.coeffs ./ p.coeffs[end])
end

Base.length(p::Polynom) = length(p.coeffs)

function Base.isless(p::Polynom{R}, q::Polynom{S}) where {R, S}
    if length(p) < length(q)
        true
    else
        false
    end
end

function (p::Polynom{R})(x::S) where {R, S}
    T = promote_type(R, S)
    res = zero(T)

    xⁱ = one(S)
    for i in 1:length(p)
        res += p.coeffs[i] * xⁱ
        xⁱ *= x
    end
    res
end

function Base.:+(p::Polynom{R}, q::Polynom{S}) where {R, S}
    T = promote_type(R, S)
    min_l = min(p, q) |> length
    
    coeffs_r = T[]
    for i in 1:min_l
        push!(coeffs_r, p.coeffs[i] + q.coeffs[i])
    end
    
    max_poly = max(p, q)
    if length(max_poly) > min_l
        for i in min_l+1:length(max_poly)
            push!(coeffs_r, max_poly.coeffs[i])
        end
    end
    Polynom{T}(coeffs_r)
end


function find_roots(p::Polynom)
    p = normalize(p)
    r₀ = r₁ = [(0.4+0.7im)^i for i in 0:length(p)-2]
    for n in 1:20
        for i in 1:length(r₀)
            r₁[i] = r₀[i] - p(r₀[i]) / prod((r₀[i] - rⱼ) for (j, rⱼ) in enumerate(r₀) if j != i)
        end
        r₀ .= r₁
        @show p.(r₀)
    end
    r₀
end