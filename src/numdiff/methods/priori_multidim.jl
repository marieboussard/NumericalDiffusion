struct PrioriMultidim{mdtype_type<:ModifiedDataType} <: QuantifMethod
    mdtype::mdtype_type
    Posteriori(mdtype::ModifiedDataType=MultiMeanMD()) = new{typeof(mdtype)}(mdtype)
end
