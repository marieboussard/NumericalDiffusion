abstract type Mesh end

struct OneDMesh <: Mesh
    Nx::Int
    xmin::Float64
    xmax::Float64 
    dx::Float64
    x::Vector{Float64}
    interfaces::Vector{Float64}

    function OneDMesh(Nx::Int, xmin::Real, xmax::Real)
        dx = (xmax - xmin) / Nx
        x = collect(LinRange(xmin+dx/2, xmax-dx/2, Nx))
        interfaces = collect(LinRange(xmin, xmax, Nx + 1))
        new(Nx, xmin, xmax, dx, x, interfaces)
    end
end

abstract type TwoDMesh <: Mesh end

struct TwoDCartesian <: TwoDMesh
    Nx::Int
    Ny::Int
    xmin::Float64
    xmax::Float64 
    ymin::Float64
    ymax::Float64 
    dx::Float64 
    dy::Float64 
    x::Vector{Float64} 
    y::Vector{Float64}
end


struct Parameters{MeshType<:Mesh}
    mesh::MeshType
    t0::Float64
    tf::Float64
    CFL_factor::Float64
end

dimension(::OneDMesh) = 1
dimension(::TwoDMesh) = 2