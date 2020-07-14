"do nothing by default"
function constraint_mc_model_voltage(pm::_PM.AbstractPowerModel, nw::Int)
end


"Generic thermal limit constraint from-side"
function constraint_mc_thermal_limit_from(pm::_PM.AbstractPowerModel, nw::Int, f_idx::Tuple{Int,Int,Int}, rate_a::Vector{<:Real})
    f_connections = ref(pm, nw, :branch, f_idx[1])["f_connections"]
    p_fr = [var(pm, nw, :p, f_idx)[c] for c in f_connections]
    q_fr = [var(pm, nw, :q, f_idx)[c] for c in f_connections]

    mu_sm_fr = JuMP.@constraint(pm.model, p_fr.^2 + q_fr.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, f_idx[1])[:mu_sm_fr] = mu_sm_fr
    end
end


"Generic thermal limit constraint to-side"
function constraint_mc_thermal_limit_to(pm::_PM.AbstractPowerModel, nw::Int, t_idx::Tuple{Int,Int,Int}, rate_a::Vector{<:Real})
    t_connections = ref(pm, nw, :branch, t_idx[1])["t_connections"]
    p_to = [var(pm, nw, :p, t_idx)[c] for c in t_connections]
    q_to = [var(pm, nw, :q, t_idx)[c] for c in t_connections]

    mu_sm_to = JuMP.@constraint(pm.model, p_to.^2 + q_to.^2 .<= rate_a.^2)

    if _IM.report_duals(pm)
        sol(pm, nw, :branch, t_idx[1])[:mu_sm_to] = mu_sm_to
    end
end


"on/off bus voltage magnitude constraint"
function constraint_mc_bus_voltage_magnitude_on_off(pm::_PM.AbstractPowerModel, n::Int, i::Int, vmin, vmax)
    vm = var(pm, n, :vm, i)
    z_voltage = var(pm, n, :z_voltage, i)

    for c in conductor_ids(pm, n)
        if isfinite(vmax[c])
            JuMP.@constraint(pm.model, vm[c] <= vmax[c]*z_voltage)
        end

        if isfinite(vmin[c])
            JuMP.@constraint(pm.model, vm[c] >= vmin[c]*z_voltage)
        end
    end
end


"on/off bus voltage magnitude squared constraint for relaxed formulations"
function constraint_mc_bus_voltage_magnitude_sqr_on_off(pm::_PM.AbstractPowerModel, n::Int, i::Int, vmin, vmax)
    w = var(pm, n, :w, i)
    z_voltage = var(pm, n, :z_voltage, i)

    for c in conductor_ids(pm, n)
        if isfinite(vmax[c])
            JuMP.@constraint(pm.model, w[c] <= vmax[c]^2*z_voltage)
        end

        if isfinite(vmin[c])
            JuMP.@constraint(pm.model, w[c] >= vmin[c]^2*z_voltage)
        end
    end
end


function constraint_mc_gen_power_setpoint_real(pm::_PM.AbstractPowerModel, n::Int, i, pg)
    pg_var = var(pm, n, :pg, i)
    JuMP.@constraint(pm.model, pg_var .== pg)
end


"on/off constraint for generators"
function constraint_mc_gen_power_on_off(pm::_PM.AbstractPowerModel, n::Int, i::Int, pmin, pmax, qmin, qmax)
    pg = var(pm, n, :pg, i)
    qg = var(pm, n, :qg, i)
    z = var(pm, n, :z_gen, i)

    for c in conductor_ids(pm, n)
        if isfinite(pmax[c])
            JuMP.@constraint(pm.model, pg[c] .<= pmax[c].*z)
        end

        if isfinite(pmin[c])
            JuMP.@constraint(pm.model, pg[c] .>= pmin[c].*z)
        end

        if isfinite(qmax[c])
            JuMP.@constraint(pm.model, qg[c] .<= qmax[c].*z)
        end

        if isfinite(qmin[c])
            JuMP.@constraint(pm.model, qg[c] .>= qmin[c].*z)
        end
    end
end


""
function constraint_mc_storage_thermal_limit(pm::_PM.AbstractPowerModel, n::Int, i, rating)
    connections = ref(pm, n, :storage, i)["connections"]
    ps = [var(pm, n, :ps, i)[c] for c in connections]
    qs = [var(pm, n, :qs, i)[c] for c in connections]

    JuMP.@constraint(pm.model, ps.^2 + qs.^2 .<= rating.^2)
end
