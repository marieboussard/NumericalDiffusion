include("../../src/include_file.jl")

Nx = 102
CFL_factor = 0.5
equation = burgers()
scheme = FVScheme(Euler(), Roe(CFL_factor))
boundsType = NormalBounds()
nb_it = 1
paramBound = 10

# counter_ex = find_counter_example(Nx, CFL_factor, equation, scheme; modifiedDataType = AsymmetricModifiedData(), boundsType=boundsType, nb_it=nb_it, paramBound=paramBound)

# @show counter_ex;

a1, b1, a2, b2 = 1.0, 2.0, 2.0, 4.5
@show counter_ex_epsilon([a1,b1,a2,b2], Nx, CFL_factor, equation, scheme, AsymmetricModifiedData(), boundsType)