abstract type QuantifParams end

struct QuantifData{dataType<:AbstractArray}
    u0::dataType
    u::dataType
end