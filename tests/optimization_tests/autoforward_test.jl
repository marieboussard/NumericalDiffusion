include("../src/include_file.jl")

function constraint_test(vh, alpha::AbstractArray{T}) where T
    #v = zeros(T,3)
    v = copy(vh)
    x = [i for i in eachindex(v)]
    for k in eachindex(v)
        v[k] += x[k] + alpha[k]
        #v[k] += alpha[k]
    end
    sum(v)
end

function cost_test(alpha::AbstractArray{T}) where T
    sum(alpha.^2)
end

function find_optimal_test()
    alphaG_init, _p = zeros(10), 0.0

    vh = zeros(10)

    cons(res, alphaG, p) = (res .= [constraint_test(vh, alphaG)])
    @show constraint_test(vh, alphaG_init)

    optprob = OptimizationFunction((x,p) -> cost_test(x), Optimization.AutoForwardDiff(); cons = cons)
    prob = OptimizationProblem(optprob, alphaG_init, _p, lcons = [-1000.0], ucons = [0.0])
    solAlphaG = solve(prob, Ipopt.Optimizer())

end

println("Test")
res = find_optimal_test()
@show "Succes!"