struct LogConfig
    ulog::Bool
    tlog::Bool
    dtlog::Bool
    fnumlog::Bool
end

const DefaultLogConfig = LogConfig(false, false, false, false)

struct LogBook
    config::LogConfig
    ulog
    tlog::Vector{Float64}
    dtlog::Vector{Float64}
    fnumlog

    function LogBook(config)
        ulog = AbstractArray{Float64}[]
        tlog = Float64[]
        dtlog = Float64[]
        fnumlog = AbstractArray{Float64}[]
        new(config, ulog, tlog, dtlog, fnumlog)
    end
end