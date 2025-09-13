# f = -1

# G(x, beta) = f^2/2*(x^3/3-x^2/2+x/4) - f*(-f*x^3/6 + f*x^2/4 + beta*x)

# xminus(beta) = -1/2 + sqrt(f/2*(f/2+4*beta))/f 
# xplus(beta) = -1/2 - sqrt(f/2*(f/2+4*beta))/f

# function J(beta)
#     if beta <= 0
#         return 0.0
#     elseif beta < -f/8
#         return G(1, beta) - G(xplus(beta), beta) + G(xminus(beta), beta)
#     else
#         return -f^2/24 - f*beta 
#     end
# end

# betavec = LinRange(-1, 1, 100)
# Jvec = J.(betavec)

f = 1

J(beta) = -f^2/24 - f*beta 
betavec = LinRange(0, 1, 100)
Jvec = J.(betavec)


using CairoMakie 

fig = Figure()
ax = Axis(fig[1,1], title="f = "*string(f), ylabel="J(beta)", xlabel="beta")
lines!(ax, betavec, Jvec)
fig