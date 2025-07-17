struct LogConfig
    ulog::Bool
    tlog::Bool
    dtlog::Bool
    fnumlog::Bool
    fcontlog::Bool
end

const DefaultLogConfig = LogConfig(false, false, false, false, false)

struct LogBook{utype, fnumtype}
    config::LogConfig
    ulog::utype
    tlog::Vector{Float64}
    dtlog::Vector{Float64}
    fnumlog::fnumtype
    fcontlog::fnumtype

    function LogBook(config, u, fnum)
        ulog = (typeof(u))[]
        tlog = (Float64)[]
        dtlog = (Float64)[]
        fnumlog = (typeof(fnum))[]
        fcontlog = (typeof(fnum))[]
        # ulog = AbstractArray{Float64}[]
        # tlog = Float64[]
        # dtlog = Float64[]
        # fnumlog = AbstractArray{Float64}[]
        #fnumlog = []
        new{typeof(ulog), typeof(fnumlog)}(config, ulog, tlog, dtlog, fnumlog, fcontlog)
    end
end