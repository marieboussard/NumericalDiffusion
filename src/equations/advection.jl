# 1D ADVECTION

struct Advection <: AbstractEquationFun 
a::Float64
end

flux(eq::Advection, u) = eq.a*u
Dflux(eq::Advection, u) = eq.a

flux!(eq::Advection, u, res) = copyto!(res, eq.a*u)
Dflux!(eq::Advection, u, res) = copyto!(res, eq.a)

u0_gauss(x::Real; xm=0.0, sigma=0.1) = 1/sqrt(2*pi*sigma)*exp(-(x-xm)^2/sigma^2)
u0_gauss(x::AbstractArray) = u0_gauss.(x)

AdvectionExample = Equation(OneD(), 1, Scalar(), Advection(2.0), u0_gauss)

# 2D ADVECTION

    # CONSTANT COEFFICIENTS

struct ConstAdvection2D <: AbstractEquationFun
    a::Float64
    b::Float64
end

flux_f(eq::ConstAdvection2D, u) = eq.a*u
flux_h(eq::ConstAdvection2D, u) = eq.b*u
Dflux_f(eq::ConstAdvection2D, u) = eq.a
Dflux_h(eq::ConstAdvection2D, u) = eq.b

# IN PLACE VERSION
# function flux_f!(eq::ConstAdvection2D, integrator)
#     integrator.fcont.fcont .= eq.a*integrator.u
# end
# function flux_h!(eq::ConstAdvection2D, integrator)
#     integrator.fcont.hcont .= eq.b*integrator.u
# end
function flux!(eq::ConstAdvection2D, u::AbstractArray, resf::AbstractArray, resh::AbstractArray)
    copyto!(resf, eq.a*u)
    copyto!(resh, eq.b*u)
end

u0_gauss2(x::Real, y::Real; xm=0.0, ym=0.0, sigmax=0.2, sigmay=0.2, A=2.0) = A*exp(-(x-xm)^2/(2*sigmax^2)-(y-ym)^2/(2*sigmay^2))

Advection2Example = Equation(TwoD(), 1, Scalar(), ConstAdvection2D(2.0, 1.0), u0_gauss2)

    # DIVERGENCE FREE VELOCITY FIELD

struct Advection2D <: AbstractEquationFun
    cnum::Array{Float64,3}
    function Advection2D(mesh::TwoDCartesian, varphi::Base.Callable=phi)
        @unpack Nx, Ny, x, y = mesh 
        cnum = zeros(Nx, Ny, 2)
        for j in 1:Nx
            for k in 1:Ny
                cnum[j,k,:] .= velfield(x[j],y[k],varphi)
            end
        end
        new(cnum)
    end
end

phi(v) = exp(-v)
velfield(x,y, varphi::Base.Callable) = varphi(x^2+y^2)*[y, -x]

flux_f(eq::Advection2D, u::Float64, j, k) = eq.cnum[j,k,1]*u
flux_h(eq::Advection2D, u::Float64, j, k) = eq.cnum[j,k,2]*u
Dflux_f(eq::Advection2D, ::Float64, j, k) = eq.cnum[j,k,1]
Dflux_h(eq::Advection2D, ::Float64, j, k) = eq.cnum[j,k,2]

function flux_f(eq::Advection2D, u::Matrix)
    a = view(eq.cnum, :, :, 1)
    a.*u
end
function flux_h(eq::Advection2D, u::Matrix)
    b = view(eq.cnum, :, :, 2)
    b.*u
end
Dflux_f(eq::Advection2D, ::Matrix) = eq.cnum[:,:,1]
Dflux_h(eq::Advection2D, ::Matrix) = eq.cnum[:,:,2]

# IN PLACE VERSION
# function flux_f!(eq::Advection2D, integrator)
#     a = view(eq.cnum, :, :, 1)
#     integrator.fcont.fcont .= a.*integrator.u
# end
# function flux_h!(eq::Advection2D, integrator)
#     b = view(eq.cnum, :, :, 2)
#     integrator.fcont.hcont .= b.*integrator.u
# end
function flux!(eq::Advection2D, u::AbstractArray, resf::AbstractArray, resh::AbstractArray)
    a = view(eq.cnum, :, :, 1)
    b = view(eq.cnum, :, :, 2)
    copyto!(resf, a.*u)
    copyto!(resh, b.*u)
end

advection2_vecfield(mesh::TwoDCartesian, varphi::Base.Callable=phi; kwargs...) = Equation(TwoD(), 1, Scalar(), Advection2D(mesh, varphi), (x,y) -> u0_gauss2(x,y;kwargs...))

# EXACT SOLUTION IN THE CASE φ(v)=1

exact_advection_sol(t, x, y, u0::Base.Callable) = u0(x*cos(t)-y*sin(t), x*sin(t)+y*cos(t))