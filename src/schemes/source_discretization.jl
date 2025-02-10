# struct SourceDiscretize{F1<:Base.Callable} <: AbstractSourceDiscretize
#     source_term::F1
# end

# struct InterpolDiscretize <: AbstractSourceDiscretize end

# source_term(source_discretize::SourceDiscretize, args...) = source_discretize.source_term(args...)

# source_term(::InterpolDiscretize, sourcefun, x) = sourcefun.(x)

abstract type SourceDiscretize{eqFunType<:AbstractEquationFun} end
struct Pointwise <: SourceDiscretize{SaintVenant} end
struct HRDisc <: SourceDiscretize{SaintVenant} end
