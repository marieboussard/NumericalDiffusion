abstract type Mesh{T<:Real} end

struct OneDMesh{T<:Real} <: Mesh{T}
    Nx::Int
    xmin::T
    xmax::T 
    dx::T
    x::Vector{T}
    interfaces::Vector{T}
end

abstract type TwoDMesh{T<:Real} <: Mesh{T} end

struct TwoDCartesian{T<:Real} <: TwoDMesh{T}
    Nx::Int
    Ny::Int
    xmin::T
    xmax::T 
    ymin::T
    ymax::T 
    dx::T 
    dy::T 
    x::Vector{T} 
    y::Vector{T}
end


struct Parameters{T<:Real, MeshType<:Mesh{T}}
    mesh::MeshType
    t0::T
    tf::T
    CFL_factor::T
end

dimension(::OneDMesh) = 1
dimension(::TwoDMesh) = 2