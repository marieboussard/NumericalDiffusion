struct LogConfig
    u::Bool
    t::Bool
    dt::Bool
end

DefaultLogConfig = LogConfig(false, false, false)

struct LogBook
    config::LogConfig
    u_log
    t_log
    dt_log

    function LogBook(config, T)
        u_log = T[]
        t_log = T[]
        dt_log = T[]
        new(config, u_log, t_log, dt_log)
    end
end