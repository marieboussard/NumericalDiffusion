# A projection algorithm

```@docs
newton_lcp(M::AbstractMatrix, q::AbstractVector, z0::AbstractVector, w0::AbstractVector; maxiter::Int=1000, c=1)
```

```@docs
compute_entropic_G(params::Parameters, equation::Equation; bound_mode=SingleBound(), weights=AbsWeights(), maxiter_uzawa=1000, maxiter_newton=1000, time_scheme::TimeScheme=Euler(), space_scheme::SpaceScheme=Rusanov(), use_newton=true)
```

```@docs
UzawaNewtonSolution{soltype<:Solution, wtype<:AbstractMatrix, atype<:AbstractMatrix, gtype<:AbstractVector, Mtype<:AbstractMatrix, ptype<:AbstractVector, bound_mode_type<:BoundMode, weights_type<:AbstractNormWeights}
```

```@docs
init_optim_components
```
