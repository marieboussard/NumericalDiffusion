# Numerical flux

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