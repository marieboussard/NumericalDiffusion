struct LogConfig
    ulog::Bool
    tlog::Bool
    dtlog::Bool
end

const DefaultLogConfig = LogConfig(false, false, false)

struct LogBook
    config::LogConfig
    ulog
    tlog
    dtlog

    function LogBook(config)
        ulog = AbstractArray{Float64}[]
        tlog = Float64[]
        dtlog = Float64[]
        new(config, ulog, tlog, dtlog)
    end
end