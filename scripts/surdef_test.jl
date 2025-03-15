u1 = zeros(6)
K1 = ones(1)

u2 = zeros(6,2)
K2 = ones(2)

# function ftest(u, K)
#     for I in CartesianIndices(u)
#         println(I)
#     end
# end

# @show ftest(u1, K1)
# @show ftest(u2, K2)

ftest(u, K) = K .= mean(u, dims=1)'
ftest2(u, K) = K .= mean(u, dims=1)

@show @allocated ftest(u1, K1)
@show @allocated ftest2(u1, K1)