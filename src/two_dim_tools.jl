# STORAGE OF u

function initialize_u(::TwoD, ::NoSource, equation::AbstractEquation, params::Parameters)
    @unpack Nx, Ny, x, y = params.mesh
    uinit = zeros(eltype(x), (Nx, Ny))
    for j in eachindex(x)
        for k in eachindex(y)
            uinit[j,k] = equation.initcond(x[j], y[k])
        end
    end
    uinit
end


# NUMERICAL FLUX

abstract type AbstractFnum{T<:EquationDim, S<:EquationType} end

struct OneDFnum{S<:EquationType, ftype<:AbstractArray} <: AbstractFnum{OneD, S}
    fnum::ftype
    function OneDFnum(equation::Equation, mesh::Mesh)
        fnum = init_fnum_one(equation.eqtype, equation, mesh)
        new{typeof(equation.eqtype), typeof(fnum)}(fnum)
    end
end

struct TwoDFnum{S<:EquationType, ftype<:AbstractArray} <: AbstractFnum{TwoD, S}
    fnum::ftype
    hnum::ftype
    function TwoDFnum(equation::Equation, mesh::Mesh)
        fnum, hnum = init_fnum_two(equation, mesh)
        new{typeof(equation.eqtype), typeof(fnum)}(fnum, hnum)
    end
end

init_fnum_one(::Scalar, equation::Equation, mesh::Mesh) = zeros(Float64, mesh.Nx+1)
init_fnum_one(::System, equation::Equation, mesh::Mesh) = zeros(Float64, (mesh.Nx+1, equation.p))

function init_fnum_two(equation::Equation, mesh::Mesh)
    fnum = zeros(Float64, mesh.Nx+1, mesh.Ny, equation.p)
    hnum = zeros(Float64, mesh.Nx, mesh.Ny+1, equation.p)
    fnum, hnum
end

# EVALUATED CONTINUOUS FLUX

abstract type AbstractFcont{T<:EquationDim, S<:EquationType} end

struct OneDFcont{S<:EquationType, ftype<:AbstractArray} <: AbstractFcont{OneD, S}
    fcont::ftype
    function OneDFcont(equation::Equation, u)
        fcont = init_fcont_one(equation, u)
        new{typeof(equation.eqtype), typeof(fcont)}(fcont)
    end
end

struct TwoDFcont{S<:EquationType, ftype<:AbstractArray} <: AbstractFcont{TwoD, S}
    fcont::ftype
    hcont::ftype
    function TwoDFcont(equation::Equation, u)
        fcont, hcont = init_fcont_two(equation, u)
        new{typeof(equation.eqtype), typeof(fcont)}(fcont, hcont)
    end
end

init_fcont_one(equation::Equation, u) = flux(equation.funcs, u)
# init_fcont_one(::System, equation::Equation, mesh::Mesh) = zeros(Float64, (mesh.Nx+1, equation.p))

function init_fcont_two(equation::Equation, u)
    fcont = flux_f(equation.funcs, u)
    hcont = flux_h(equation.funcs, u)
    fcont, hcont
end