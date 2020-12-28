using Ipopt; using PowerModels
using SortingLab

solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

network_name = "data/matpower/case14.m"
#network_name = "data/matpower/pglib_opf_case39_epri.m"
#network_name = "data/matpower/pglib_opf_case57_ieee.m"
#network_name = "data/matpower/pglib_opf_case73_ieee_rts.m"
#network_name = "data/matpower/pglib_opf_case89_pegase.m"
#network_name = "data/matpower/pglib_opf_case118_ieee.m"
#network_name = "data/matpower/pglib_opf_case162_ieee_dtc.m"
#network_name = "data/matpower/pglib_opf_case179_goc.m"
#network_name = "data/matpower/pglib_opf_case200_activ.m"   
#network_name = "data/matpower/pglib_opf_case240_pserc.m" 

network_data = PowerModels.parse_file(network_name)
network_data_new = PowerModels.parse_file(network_name)
network_data_new_v1 = PowerModels.parse_file(network_name)
network_data_new_v2 = PowerModels.parse_file(network_name)
network_data_new_v3 = PowerModels.parse_file(network_name)
network_data_new_v4 = PowerModels.parse_file(network_name)
network_data_new_v5 = PowerModels.parse_file(network_name)
network_data_new_v6 = PowerModels.parse_file(network_name)

s = Dict("output" => Dict("branch_flows" => true) ) 

result_ac_original = run_ac_opf(network_data, with_optimizer(Ipopt.Optimizer); setting = s)

# Check for number of generators in the network
g_upper = length(network_data_new["gen"])

# for i=1:g_upper

#     g=string(i)
#     network_data["gen"][g]["pg"] = result_soc["solution"]["gen"][g]["pg"]
#     network_data["gen"][g]["qg"] = result_soc["solution"]["gen"][g]["qg"]

# end

bus_id = [bus["index"] for (i,bus) in network_data_new["gen"]]

Gen_active_orig = [bus["pg"] for (i,bus) in network_data_new["gen"]]
Gen_reactive_orig = [bus["qg"] for (i,bus) in network_data_new["gen"]]

p1 = plot(bus_id, Gen_active_orig, seriestype = :scatter) # Make a line plot
p3 = plot(bus_id,Gen_reactive_orig , seriestype = :scatter , xlabel = "reactive power", lw = 3, title = "modified")

for (i,bus) in result_ac_original["solution"]["gen"]
    network_data_new["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
    network_data_new["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])

    network_data_new_v1["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
    network_data_new_v1["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])

    network_data_new_v2["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
    network_data_new_v2["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])

    network_data_new_v3["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
    network_data_new_v3["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])

    network_data_new_v4["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
    network_data_new_v4["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])

    network_data_new_v5["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
    network_data_new_v5["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])

    network_data_new_v6["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
    network_data_new_v6["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])

end

Gen_active_new = [bus["pg"] for (i,bus) in network_data_new["gen"]]
Gen_reactive_new = [bus["qg"] for (i,bus) in network_data_new["gen"]]


p2 = plot!(p1, bus_id,Gen_active_new, seriestype = :scatter) # Make a scatter plot
p4 = plot!(p3, bus_id, Gen_reactive_new, seriestype = :scatter) # Four histograms each with 10 points? Why not!
plot(p2, p4, layout = (1, 2), legend = false)

load_bus_id = [bus["index"] for (i,bus) in network_data_new["load"]]



for (i,bus) in network_data["load"]
    network_data_new["load"][i]["pd"] = 1.25*deepcopy( network_data["load"][i]["pd"] )
    network_data_new["load"][i]["qd"] = 1.25*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v1["load"][i]["pd"] = 1.2*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v1["load"][i]["qd"] = 1.2*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v2["load"][i]["pd"] = 1.1*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v2["load"][i]["qd"] = 1.1*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v3["load"][i]["pd"] = 1.05*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v3["load"][i]["qd"] = 1.05*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v4["load"][i]["pd"] = 0.95*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v4["load"][i]["qd"] = 0.95*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v5["load"][i]["pd"] = 0.9*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v5["load"][i]["qd"] = 0.9deepcopy(network_data["load"][i]["qd"])

    network_data_new_v6["load"][i]["pd"] = 0.8*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v6["load"][i]["qd"] = 0.8*deepcopy(network_data["load"][i]["qd"])
end
# Active_orig = [bus["pd"] for (i,bus) in network_data["load"]]
# Reactive_orig = [bus["qd"] for (i,bus) in network_data["load"]]

# Active_new = [bus["pd"] for (i,bus) in network_data_new["load"]]
# Reactive_new = [bus["qd"] for (i,bus) in network_data_new["load"]]

# Active_new_v1 = [bus["pd"] for (i,bus) in network_data_new_v1["load"]]
# Reactive_new_v1 = [bus["qd"] for (i,bus) in network_data_new_v1["load"]]

# Active_new_v2 = [bus["pd"] for (i,bus) in network_data_new_v2["load"]]
# Reactive_new_v2 = [bus["qd"] for (i,bus) in network_data_new_v2["load"]]

# Active_new_v3 = [bus["pd"] for (i,bus) in network_data_new_v3["load"]]
# Reactive_new_v3 = [bus["qd"] for (i,bus) in network_data_new_v3["load"]]

# Active_new_v4 = [bus["pd"] for (i,bus) in network_data_new_v4["load"]]
# Reactive_new_v4 = [bus["qd"] for (i,bus) in network_data_new_v4["load"]]

# Active_new_v5 = [bus["pd"] for (i,bus) in network_data_new_v5["load"]]
# Reactive_new_v5 = [bus["qd"] for (i,bus) in network_data_new_v5["load"]]

# Active_new_v6 = [bus["pd"] for (i,bus) in network_data_new_v6["load"]]
# Reactive_new_v6 = [bus["qd"] for (i,bus) in network_data_new_v6["load"]]

# dict_active_power= [load_bus_id Active_orig Active_new Active_new_v1 Active_new_v2 Active_new_v3 Active_new_v4 Active_new_v5 Active_new_v6]
# dict_REactive_power= [load_bus_id Reactive_orig Reactive_new Reactive_new_v1 Reactive_new_v2 Reactive_new_v3 Reactive_new_v4 Reactive_new_v5 Reactive_new_v6]

# idx_acti = fsortperm(dict_active_power[:,1])
# Act_sorted = dict_active_power[idx_acti,:]

# idx_REacti = fsortperm(dict_REactive_power[:,1])
# React_sorted = dict_REactive_power[idx_REacti,:]

# x_ac = Act_sorted[:,1]
# y_ac = Act_sorted[:,2:end]
# y_re = React_sorted[:,2:end]

# s1 = plot(x_ac, y_ac, seriestype = :scatter,  title = "active load") # Make a line plot
# s2 = plot(x_ac,y_re , seriestype = :scatter , xlabel = "bus id", lw = 3, title = "reactive load")
# plot(s1, s2, layout = (1, 2), legend = false)



l0 = PowerModels.component_table(network_data,"load", ["pd", "qd"])
l1 = PowerModels.component_table(network_data_new,"load", ["pd", "qd"])
l2 = PowerModels.component_table(network_data_new_v1,"load", ["pd", "qd"])
l3 = PowerModels.component_table(network_data_new_v2,"load", ["pd", "qd"])
l4 = PowerModels.component_table(network_data_new_v3,"load", ["pd", "qd"])
l5 = PowerModels.component_table(network_data_new_v4,"load", ["pd", "qd"])
l6 = PowerModels.component_table(network_data_new_v5,"load", ["pd", "qd"])
l7 = PowerModels.component_table(network_data_new_v6,"load", ["pd", "qd"])

load_active_matrix = Array{Float64}([l0[:,1] l0[:,2] l1[:,2] l2[:,2] l3[:,2] l4[:,2] l5[:,2] l6[:,2] l7[:,2]])
load_REactive_matrix = Array{Float64}([l0[:,1] l0[:,3] l1[:,3] l2[:,3] l3[:,3] l4[:,3] l5[:,3] l6[:,3] l7[:,3]])
x_load= load_active_matrix[:,1]
y_active= load_active_matrix[:,2:end]
y_REactive= load_REactive_matrix[:,2:end]

plot(x_load, y_active)
plot(x_load, y_REactive)

active_matrix_difference = Array{Float64}([l1[:,1] l0[:,2]-l2[:,2] l0[:,2]-l3[:,2] l0[:,2]-l4[:,2] l0[:,2]-l5[:,2] l0[:,2]-l6[:,2] l0[:,2]-l7[:,2]])
reactive_difference = Array{Float64}([l1[:,1] l0[:,3]-l2[:,3] l0[:,3]-l3[:,3] l0[:,3]-l4[:,3] l0[:,3]-l5[:,3] l0[:,3]-l6[:,3] l0[:,3]-l7[:,3]])

plot(active_matrix_difference[:,1], active_matrix_difference[:,2:end])

plot(active_matrix_difference[:,1], reactive_difference[:,2:end])

# s1 = plot(load_bus_id, Active_orig, seriestype = :scatter,  title = "active load") # Make a line plot
# s2 = plot(load_bus_id,Reactive_orig , seriestype = :scatter , xlabel = "bus id", lw = 3, title = "reactive load")

# s3 = plot!(s1,load_bus_id, Active_new, seriestype = :scatter)
# s4 = plot!(s2,load_bus_id,Reactive_new , seriestype = :scatter )

# s5 = plot!(s3,load_bus_id, Active_new_v1, seriestype = :scatter)
# s6 = plot!(s4,load_bus_id,Reactive_new_v1 , seriestype = :scatter )

# s7 = plot!(s5,load_bus_id, Active_new_v1, seriestype = :scatter)
# s8 = plot!(s6,load_bus_id,Reactive_new_v1 , seriestype = :scatter )

# s9 = plot!(s7,load_bus_id, Active_new_v1, seriestype = :scatter)
# s10 = plot!(s8,load_bus_id,Reactive_new_v1 , seriestype = :scatter )

# s11 = plot!(s9,load_bus_id, Active_new_v1, seriestype = :scatter)
# s12 = plot!(s10,load_bus_id,Reactive_new_v1 , seriestype = :scatter )

# s13 = plot!(s11,load_bus_id, Active_new_v1, seriestype = :scatter)
# s14 = plot!(s12,load_bus_id,Reactive_new_v1 , seriestype = :scatter )

# s15 = plot!(s13,load_bus_id, Active_new_v1, seriestype = :scatter)
# s16 = plot!(s14,load_bus_id,Reactive_new_v1 , seriestype = :scatter )

# plot(s15, s16, layout = (1, 2), legend = false)


result_new = run_ac_opf(network_data_new, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v1 = run_ac_pf(network_data_new_v1, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v2 = run_ac_pf(network_data_new_v2, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v3 = run_ac_pf(network_data_new_v3, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v4 = run_ac_pf(network_data_new_v4, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v5 = run_ac_pf(network_data_new_v5, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v6 = run_ac_pf(network_data_new_v6, with_optimizer(Ipopt.Optimizer); setting = s)


# V_ac_orig = [bus["vm"] for (i,bus) in result_ac_original["solution"]["bus"]]
# V_ac_new= [bus["vm"] for (i,bus) in result_new["solution"]["bus"]]
# V_ac_new_v1= [bus["vm"] for (i,bus) in result_new_v1["solution"]["bus"]]
# V_ac_new_v2= [bus["vm"] for (i,bus) in result_new_v2["solution"]["bus"]]
# V_ac_new_v3= [bus["vm"] for (i,bus) in result_new_v3["solution"]["bus"]]
# V_ac_new_v4= [bus["vm"] for (i,bus) in result_new_v4["solution"]["bus"]]
# V_ac_new_v5= [bus["vm"] for (i,bus) in result_new_v5["solution"]["bus"]]
# V_ac_new_v6= [bus["vm"] for (i,bus) in result_new_v6["solution"]["bus"]]

# V_ac_orig_ang = [bus["va"] for (i,bus) in result_ac_original["solution"]["bus"]]
# V_ac_new_ang= [bus["va"] for (i,bus) in result_new["solution"]["bus"]]
# V_ac_new_v1_ang= [bus["va"] for (i,bus) in result_new_v1["solution"]["bus"]]
# V_ac_new_v2_ang= [bus["va"] for (i,bus) in result_new_v2["solution"]["bus"]]
# V_ac_new_v3_ang= [bus["va"] for (i,bus) in result_new_v3["solution"]["bus"]]
# V_ac_new_v4_ang= [bus["va"] for (i,bus) in result_new_v4["solution"]["bus"]]
# V_ac_new_v5_ang= [bus["va"] for (i,bus) in result_new_v5["solution"]["bus"]]
# V_ac_new_v6_ang= [bus["va"] for (i,bus) in result_new_v6["solution"]["bus"]]

# V_upper_lim = [bus["vmax"] for (i,bus) in network_data["bus"]]
# V_lower_lim = [bus["vmin"] for (i,bus) in network_data["bus"]]

# plot(V_lower_lim, label = "Lower limit of voltage")
# plot!(V_upper_lim, label = "Upper limit of voltage")
# plot!(V_ac_orig, label = "AC power flow")
# plot!(V_ac_new, label = "AC-PF original")
# plot!(V_ac_new_v1, label = "v1")
# plot!(V_ac_new_v2, label = "v2")
# plot!(V_ac_new_v3, label = "v3")
# plot!(V_ac_new_v4, label = "v4")
# plot!(V_ac_new_v5, label = "v5")
# plot!(V_ac_new_v6, label = "v6")

# plot(V_ac_orig_ang, label = "AC power flow")
# plot!(V_ac_new_ang, label = "AC-PF original")
# plot!(V_ac_new_v1_ang, label = "v1")
# plot!(V_ac_new_v2_ang, label = "v2")
# plot!(V_ac_new_v3_ang, label = "v3")
# plot!(V_ac_new_v4_ang, label = "v4")
# plot!(V_ac_new_v5_ang, label = "v5")
# plot!(V_ac_new_v6_ang, label = "v6")

thermal_limit= ( [branch["rate_a"] for (i,branch) in network_data["branch"]] )
#thermal_ac_orig  = ( abs.( [branch["pt"] for (i,branch) in result_ac_original["solution"]["branch"]] ) )

#loss_ac =  Dict(name => data["p_to"]+data["p_from"] for (name, data) in result_ac_original["solution"]["branch"])
#loss_dc =  Dict(name => data["p_to"]+data["p_from"] for (name, data) in result_ac_original["solution"]["dcline"])

#p_res = PowerModels.print_summary(result_ac_original["solution"])


#Dict(name => data["va"] for (name, data) in result_ac_original["solution"]["bus"])








m1 = PowerModels.component_table(network_data, "bus", ["vmin", "vmax"])
m0 = PowerModels.component_table(result_ac_original["solution"], "bus", ["vm", "va"])
m01 = PowerModels.component_table(result_new["solution"], "bus", ["vm", "va"])
m2 = PowerModels.component_table(result_new_v1["solution"], "bus", ["vm", "va"])
m3 = PowerModels.component_table(result_new_v2["solution"], "bus", ["vm", "va"])
m4 = PowerModels.component_table(result_new_v3["solution"], "bus", ["vm", "va"])
m5 = PowerModels.component_table(result_new_v4["solution"], "bus", ["vm", "va"])
m6 = PowerModels.component_table(result_new_v5["solution"], "bus", ["vm", "va"])
m7 = PowerModels.component_table(result_new_v6["solution"], "bus", ["vm", "va"])

m8 = PowerModels.component_table(network_data, "branch", ["angmin", "angmax"])

voltage_matrix = Array{Float64}([m1 m0[:,2] m01[:,2] m2[:,2] m3[:,2] m4[:,2] m5[:,2] m6[:,2] m7[:,2]])
voltage_angle = Array{Float64}([m1[:,1] m0[:,3] m01[:,3] m2[:,3] m3[:,3] m4[:,3] m5[:,3] m6[:,3] m7[:,3]])
x_volt= voltage_matrix[:,1]
y_volt= voltage_matrix[:,2:end]
plot(x_volt, y_volt)

voltage_matrix_difference = Array{Float64}([m1[:,1] m0[:,2]-m01[:,2] m0[:,2]-m2[:,2] m0[:,2]-m3[:,2] m0[:,2]-m4[:,2] m0[:,2]-m5[:,2] m0[:,2]-m6[:,2] m0[:,2]-m7[:,2]])
voltage_angle_difference = Array{Float64}([m1[:,1] m0[:,3]-m01[:,3] m0[:,3]-m2[:,3] m0[:,3]-m3[:,3] m0[:,3]-m4[:,3] m0[:,3]-m5[:,3] m0[:,3]-m6[:,3] m0[:,3]-m7[:,3]])

plot(voltage_matrix_difference[:,1], voltage_matrix_difference[:,2:end])

plot(voltage_angle_difference[:,1], voltage_angle_difference[:,2:end])



# Bound tightening strategy



PowerModels.update_data!(network_data_new_v1, result_new_v1["solution"])
dc_flows_v1 = PowerModels.calc_branch_flow_dc(network_data_new_v1)
ac_flows_v1 = PowerModels.calc_branch_flow_ac(network_data_new_v1)
ac_flows_v2 = PowerModels.calc_branch_flow_ac(network_data_new_v2)
ac_flows_v3 = PowerModels.calc_branch_flow_ac(network_data_new_v3)
ac_flows_v4 = PowerModels.calc_branch_flow_ac(network_data_new_v4)
ac_flows_v5 = PowerModels.calc_branch_flow_ac(network_data_new_v5)
ac_flows_v6 = PowerModels.calc_branch_flow_ac(network_data_new_v6)
ac_flows_limit = PowerModels.calc_branch_flow_ac(network_data_new)

branc_details = PowerModels.component_table(network_data, "branch", ["f_bus", "t_bus", "rate_a","br_status","br_r", "br_x"])

Injection_ac_orig = PowerModels.component_table(result_ac_original["solution"],"branch",["qf", "pf", "qt","pt"])
Injection_ac_limit = PowerModels.component_table(ac_flows_limit,"branch",["qf", "pf", "qt","pt"])
Injection_ac_v1 = PowerModels.component_table(ac_flows_v1,"branch",["qf", "pf", "qt","pt"])
Injection_ac_v2 = PowerModels.component_table(ac_flows_v2,"branch",["qf", "pf", "qt","pt"])
Injection_ac_v3 = PowerModels.component_table(ac_flows_v3,"branch",["qf", "pf", "qt","pt"])
Injection_ac_v4 = PowerModels.component_table(ac_flows_v4,"branch",["qf", "pf", "qt","pt"])
Injection_ac_v5 = PowerModels.component_table(ac_flows_v5,"branch",["qf", "pf", "qt","pt"])
Injection_ac_v6 = PowerModels.component_table(ac_flows_v6,"branch",["qf", "pf", "qt","pt"])


thermal_ac_original = max( sqrt.(Injection_ac_orig[:,2].^2 + Injection_ac_orig[:,3].^2 ), sqrt.(Injection_ac_orig[:,4].^2 + Injection_ac_orig[:,5].^2 ))
thermal_ac_limit = max( sqrt.(Injection_ac_limit[:,2].^2 + Injection_ac_limit[:,3].^2 ), sqrt.(Injection_ac_limit[:,4].^2 + Injection_ac_limit[:,5].^2 ))
thermal_ac_v1 = max( sqrt.(Injection_ac_v1[:,2].^2 + Injection_ac_v1[:,3].^2 ), sqrt.(Injection_ac_v1[:,4].^2 + Injection_ac_v1[:,5].^2 ))
thermal_ac_v2 = max( sqrt.(Injection_ac_v2[:,2].^2 + Injection_ac_v2[:,3].^2 ), sqrt.(Injection_ac_v2[:,4].^2 + Injection_ac_v2[:,5].^2 ))
thermal_ac_v3 = max( sqrt.(Injection_ac_v3[:,2].^2 + Injection_ac_v3[:,3].^2 ), sqrt.(Injection_ac_v3[:,4].^2 + Injection_ac_v3[:,5].^2 ))
thermal_ac_v4 = max( sqrt.(Injection_ac_v4[:,2].^2 + Injection_ac_v4[:,3].^2 ), sqrt.(Injection_ac_v4[:,4].^2 + Injection_ac_v4[:,5].^2 ))
thermal_ac_v5 = max( sqrt.(Injection_ac_v5[:,2].^2 + Injection_ac_v5[:,3].^2 ), sqrt.(Injection_ac_v5[:,4].^2 + Injection_ac_v5[:,5].^2 ))
thermal_ac_v6 = max( sqrt.(Injection_ac_v6[:,2].^2 + Injection_ac_v6[:,3].^2 ), sqrt.(Injection_ac_v6[:,4].^2 + Injection_ac_v6[:,5].^2 ))

losses_ac_orig = 100*abs.( Injection_ac_orig[:,3] + Injection_ac_orig[:,5] )./branc_details[:,4]
losses_ac_limit = 100*abs.( Injection_ac_limit[:,3] + Injection_ac_limit[:,5] )./branc_details[:,4]
losses_ac_v1 = 100*abs.( Injection_ac_v1[:,3] + Injection_ac_v1[:,5] )./branc_details[:,4]
losses_ac_v2 = 100*abs.( Injection_ac_v2[:,3] + Injection_ac_v2[:,5] )./branc_details[:,4]
losses_ac_v3 = 100*abs.( Injection_ac_v3[:,3] + Injection_ac_v3[:,5] )./branc_details[:,4]
losses_ac_v4 = 100*abs.( Injection_ac_v4[:,3] + Injection_ac_v4[:,5] )./branc_details[:,4]
losses_ac_v5 = 100*abs.( Injection_ac_v5[:,3] + Injection_ac_v5[:,5] )./branc_details[:,4]
losses_ac_v6 = 100*abs.( Injection_ac_v6[:,3] + Injection_ac_v6[:,5] )./branc_details[:,4]

thr_loading_ac_orig = 100*thermal_ac_original./branc_details[:,4]
thr_loading_limit = 100*thermal_ac_limit./branc_details[:,4]
thr_loading_v1 = 100*thermal_ac_v1./branc_details[:,4]
thr_loading_v2 = 100*thermal_ac_v2./branc_details[:,4]
thr_loading_v3 = 100*thermal_ac_v3./branc_details[:,4]
thr_loading_v4 = 100*thermal_ac_v4./branc_details[:,4]
thr_loading_v5 = 100*thermal_ac_v5./branc_details[:,4]
thr_loading_v6 = 100*thermal_ac_v6./branc_details[:,4]

x_thermal = Injection_ac_orig[:,1]
y_thermal = [thr_loading_ac_orig thr_loading_limit thr_loading_v1 thr_loading_v2 thr_loading_v3 thr_loading_v4 thr_loading_v5 thr_loading_v6]
y_loss = [losses_ac_orig losses_ac_limit losses_ac_v1 losses_ac_v2 losses_ac_v3 losses_ac_v4 losses_ac_v5 losses_ac_v6]

plot(x_thermal, y_thermal)

plot(x_thermal, y_loss)



Generator_limits = PowerModels.component_table(network_data,"gen",["gen_bus", "pmax", "pmin", "qmax","pmin"])
Generator_output = PowerModels.component_table(result_ac_original["solution"],"gen",["pg", "qg"])
Generator_output_infla = PowerModels.component_table(result_new["solution"],"gen",["pg", "qg"])
Generator_output_v1= PowerModels.component_table(result_new_v1["solution"],"gen",["pg", "qg"])
Generator_output_v6= PowerModels.component_table(result_new_v6["solution"],"gen",["pg", "qg"])


x_gen = Generator_limits[:,1]
y_gen_act = [Generator_limits[:,3] Generator_limits[:,4] Generator_output[:,2] Generator_output_infla[:,2] Generator_output_v1[:,2] Generator_output_v6[:,2]]

y_gen_REact = [Generator_limits[:,5] Generator_limits[:,6] Generator_output[:,3] Generator_output_infla[:,3] Generator_output_v1[:,3] Generator_output_v6[:,3]]

plot(x_gen, y_gen_act)

plot(x_gen, y_gen_REact)
