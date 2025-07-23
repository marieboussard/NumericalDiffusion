using NumericalDiffusion

Nx = 20
xmin, xmax = [0, 1]
mesh = OneDMesh(Nx, xmin, xmax)

t0, tf = 0.0, 0.2
CFL_factor = 0.5
params = Parameters(mesh, t0, tf, CFL_factor)

# uinit = [1.0, 1.0, 1.0, 2.0, 2.0]

u0(x::Real) = x>0.7 ? 2.0 : 1.0
u0(x::AbstractArray) = u0.(x)

@show uinit = u0(mesh.x)

equation = Equation(OneD(), 1, Scalar(), Burgers(), u0)

@show dt = CFL_factor*mesh.dx/maximum(abs.(uinit))


space_scheme = Rusanov()
rcache = RusanovCache(uinit)

# time_scheme = Euler()
# time_cache = EulerCache(rcache)

time_scheme = RK2()
time_cache = RK2Cache(space_scheme, uinit, OneD(), params, dt)


fnum = zero(uinit)

#global_numflux!(time_scheme, space_scheme, time_cache, equation, uinit, fnum);



### MODIFIED DATA #####
mdtype = AsymmetricMD()
bounds = DefaultBounds()
dim = OneD()
eqtype = Scalar()
sL, sR = get_sL(time_scheme, space_scheme), get_sR(time_scheme, space_scheme)

indices = init_indices(mdtype, bounds, sL, sR)
utilde = init_utilde(mdtype, bounds, dim, eqtype, uinit, sL, sR)
uhat = init_uhat(mdtype, bounds, dim, eqtype, uinit, sL, sR)
ftilde = init_ftilde(mdtype, bounds, dim, eqtype, uinit, sL, sR)

time_cache = init_cache(time_scheme, space_scheme, utilde, equation.dim, params, dt)

j = 14
#utilde!()
indices .= mod1.(j-2*sL-sR+2 : j+sL+2*sR-1, Nx)
@views ushort = selectdim(uinit, 1, indices)
if ndims(uinit)==1 # Scalar equations
    for i in 1:sL+sR-1
        utilde[i] = ushort[sL+sR]
    end
    for i in sL+sR:2*(sL+sR)-1
        utilde[i] = ushort[i]
    end
    for i in 2*(sL+sR):3*(sL+sR)-2
        utilde[i] = ushort[2*(sL+sR)-1]
    end
end

@show utilde


global_numflux!(time_scheme, space_scheme, time_cache, equation, utilde, ftilde, sL, 3*sL+2*sR-2, params.mesh.Nx, -sL+1)

@show ftilde

for i in 1:2*(sL+sR)-2
    for r in 1:equation.p
        uhat[i,r] = utilde[i+sL,r] - dt/mesh.dx * (ftilde[i+1,r] - ftilde[i,r])
        # uhat[i,r] = utilde[i+sL_space,r] - dt/dx * (ftilde[i+1,r] - ftilde[i,r])
    end
end

@show uhat