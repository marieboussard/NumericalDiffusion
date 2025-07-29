# Quantifying Diffusion

```@docs
quantify_diffusion(sol::Solution, method::QuantifMethod, i=0; name="", kwargs...)
```

There are two main types of methods: A Posteriori and a Priori.
In a posteriori methods, we compute bounds on the numerical entropy flux, whereas a priori methods compute directly bounds on the numerical diffusion.

```@docs
compute_G_bounds!(estimator::Estimator)
```