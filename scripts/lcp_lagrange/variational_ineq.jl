using LinearAlgebra

function newton_lcp(Q::AbstractMatrix, Id::AbstractMatrix, b::AbstractVector, p0::AbstractVector, y0::AbstractVector; niter_max::Int=10, c=1)

    # Initiate stopping criteria for the loop
    test = -1 # Striclty negative if In != Inprev, 0 else
    niter = 0

    # Initialize variables and indices
    p = copy(p0)
    y = copy(y0)
    In = y .- c .*p .<= 0
    Ac = .!In

    while niter < niter_max && (minimum(p) < -1e-16 || minimum(y) < -1e-16)#  maximum(Q*p .- y .- b) > 1e-16    #test<0 && 

        # First update the zero part
        p[Ac] .= 0
        y[In] .= 0

        # Solving a system for the remaining indices
        Mat = [Q[:, In] -Id[:, Ac]]
        X = inv(Mat) * b
        #@show X

        # @show p[In]
        # @show X[1:sum(In)]

        p[In] .= X[1:sum(In)]
        y[Ac] .= X[sum(In) + 1 : end]

        Inprev = In

        # Update indices
        @. In = y - c *p <= 0
        @. Ac = .!In

        # Testing stopping criterion
        if sum(In) == sum(Inprev)
            test = - norm((In .-Inprev) .> 0)
        else
            test = -1
        end

        niter += 1

    end

    p, y, niter

end

# Problem definition
m=10
#Q0 = rand(m,m)

# Matrix for the upper bound only
Q0 = zeros(m,m)
for j in 1:m
    Q0[j,j] = 1
    for i in 1:m
        if j == i+1
            Q0[i,j] = -1
        end
    end    
end
Q0[m,1] = -1


Q = Q0'*Q0
Id = Matrix{eltype(Q)}(I,m,m)
b = 1 .- 2 .*rand(m)
#b = -rand(m)
c=1

# Initialization of variables
eps = 0.1
p0 = inv(Q .+ eps.*Id)*b
y0 = zero(p0)
# p0 = zero(b)
# y0 = -b

# Other parameters
niter_max = 20

# Solving with Newton method
p, y, niter = newton_lcp(Q, Id, b, p0, y0; niter_max=niter_max, c=c)

@show Q*p .- y .- b
@show p
@show y 
@show p' * y