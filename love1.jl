using Plots
using Roots

ζ(ω::Real, k::Real, ρ::Real, μ::Real) = sqrt(Complex(ρ * ω^2 / μ - k^2))

function detlove1(ω::Real, k::Real, μ1::Real, ρ1::Real, μ2::Real, ρ2::Real, H::Real)
    ζ1 = ζ(ω, k, ρ1, μ1)
    ζ2 = ζ(ω, k, ρ2, μ2)
    #println("ζ1 ", ζ1)
    #println("ζ2 ", ζ2)
    det = 1im * μ2 * ζ2 * cos(ζ1 * H) + 
        μ1 * ζ1 * sin(ζ1 * H)
    return det
end

function test()
    ω = 0.1
    ρ1 = 3000
    v1 = 3000 
    μ1 = ρ1 * v1^2
    
    ρ2 = 3500
    v2 = 3500
    μ2 = ρ2 * v2^2

    H = 30000

    c = collect(range(2000, 4000, step=10.0))
    k = ω ./ c
    
    det = [detlove1(ω, k1, μ1, ρ1, μ2, ρ2, H) for k1 in k]
    #println(k)
    #println(c)
    #println(real(det))
    #print(typeof(c), typeof(real(det)))
    
    global i_bk = 0
    for i in 1:(length(det)-1)
        if real(det[i]) * real(det[i+1]) < 0
            println("[", c[i], ", ", c[i+1], "]")
            global i_bk = i
            break
        end
    end

    f(c) = real(detlove1(ω, ω/c, μ1, ρ1, μ2, ρ2, H))

    cLove = find_zero(f, (c[i_bk], c[i_bk+1]), Bisection())
    println(cLove)
    gr()
    plot(c, real.(det),
        color = :black,
        xaxis = ("Phase velocity (m/s)"),
        yaxis = ("Determinant"),
        title = ("Searching for Love wave phase velocity at $(ω) Hz"),
        lab="Determinant",
        legend = :bottomright,
        )
    
    scatter!([cLove], [0.0], m=(:red, 0.8, Plots.stroke(1,:black)), ms=5,
        lab="Roots - Phase velocity")

    savefig("LovePhase.png")
end

test()