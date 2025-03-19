u1 = zeros(6)
K1 = ones(1)

u2 = zeros(6,2)
K2 = ones(1,2)

# function ftest(u, K)
#     for I in CartesianIndices(u)
#         println(I)
#     end
# end

# @show ftest(u1, K1)
# @show ftest(u2, K2)

# ftest(u, K) = K .= mean(u, dims=1)'
# ftest2(u, K) = K. = mean(u, dims=1)

computeK!(::MeanMD, u::AbstractMatrix, res::AbstractVector) = res .= vec(mean(u, dims=1))
computeK!(::MeanMD, u::AbstractMatrix, res::AbstractMatrix) = res .= mean(u, dims=1)

@show @allocated computeK!(MeanMD(), u2, K2)

# for i in 1:size(u2)[1]
#     @show @allocated u2[i,1] = K2[1,1]
#     @show @allocated u2[i,1] = 3.0
# end
