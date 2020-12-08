using Ipopt; using PowerModels

solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

network_name = "data/matpower/case24.m"

network_data = PowerModels.parse_file(network_name)
network_d_mod = network_data
s = Dict("output" => Dict("branch_flows" => true) ) 

case, stats = run_obbt_opf!(network_name, solver)

V_upper_lim = [bus["vmax"] for (i,bus) in network_data["bus"]]
V_lower_lim = [bus["vmin"] for (i,bus) in network_data["bus"]]

V_upper_lim_bt = [bus["vmax"] for (i,bus) in case["bus"]]
V_lower_lim_bt = [bus["vmin"] for (i,bus) in case["bus"]]

#[bus["vmax"] for (i,bus) in network_d_mod["bus"]] = V_upper_lim_bt
#[bus["vmin"] for (i,bus) in network_d_mod["bus"]] = V_lower_lim_bt

for (i,bus) in case["bus"]
    network_d_mod["bus"][i]["vmax"] = case["bus"][i]["vmax"]
    network_d_mod["bus"][i]["vmin"] = case["bus"][i]["vmin"]
end


Active_orig = [bus["pd"] for (i,bus) in network_data["load"]]
Reactive_orig = [bus["qd"] for (i,bus) in network_data["load"]]

for (i,bus) in network_d_mod["load"]
    network_d_mod["load"][i]["pd"] = 1.1*network_data["load"][i]["pd"]
    network_d_mod["load"][i]["qd"] = 1.1*network_data["load"][i]["qd"]
end

Active_new = [bus["pd"] for (i,bus) in network_data["load"]]
Reactive_new = [bus["qd"] for (i,bus) in network_data["load"]]

stats

#open("data/matpower/pglib_opf_case5_pjm-bt.m", "w") do io
#    export_matpower(io, case)
#end

result_ac_original = run_ac_opf(network_data, with_optimizer(Ipopt.Optimizer); setting = s)

result_ac_bt = run_ac_opf(case, with_optimizer(Ipopt.Optimizer); setting = s)

result_ac_bt2 = run_ac_opf(network_d_mod, with_optimizer(Ipopt.Optimizer); setting = s)

V_orig = [bus["vm"] for (i,bus) in result_ac_original["solution"]["bus"]]
V_bt = [bus["vm"] for (i,bus) in result_ac_bt["solution"]["bus"]]

V_case = [bus["vm"] for (i,bus) in case["bus"]]

plot(V_lower_lim, label = "Lower limit of voltage")
plot!(V_upper_lim, label = "Upper limit of voltage")
plot!(V_orig, label = "AC OPF original")
plot!(V_upper_lim_bt, label = "New Upper limi",  lw = 3)
plot!(V_lower_lim_bt, label = "New Lower Limit", lw = 3)
plot!(V_bt, label = "bound tightened")
plot!(V_case, label = "OBBT V output")

plot(V_upper_lim_bt, label = "New Upper limi",  lw = 3)
plot!(V_lower_lim_bt, label = "New Lower Limit", lw = 3)
plot!(V_case, label = "OBBT V output")


the_limit= sort( [branch["rate_a"] for (i,branch) in network_data["branch"]] )

the_case_bt= sort( [branch["rate_a"] for (i,branch) in case["branch"]] )

the_ac = sort(abs.( [branch["pt"] for (i,branch) in result_ac_original["solution"]["branch"]])  )

the_bt  = sort( abs.( [branch["pt"] for (i,branch) in result_ac_bt["solution"]["branch"]] ) )

plot(the_limit, label = "Matpower limit")
plot!(the_case_bt, label = "Case BT limit")
plot!(the_ac, label = "AC OPF original")
#plot!(the_bt, label = "AC OPF BT")


a1 = stats["initial_relaxation_objective"]
a2 = stats["final_relaxation_objective"]

a3 = result_ac_original["objective"]
a4 = result_ac_bt["objective"]

r_mat = [a1 a2 a3 a4]

