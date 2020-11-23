using PowerModels
using Ipopt
using Plots

# PowerFlow with PowerModels in Julia Modes 
# ACPPowerModel- AC in polar coordinates 
# ACRPowerModel- AC in rectangular coordinates 
# ACTPowerModel- AC in the w-theta space
# DCPPowerModel- DC approximation 
# DCPLLPowerModel- DC approx. with line losses 
# SOCWRPowerModel- SOC relaxation 
# QCWRPowerModel- QC relaxation

network_name = "data/matpower/case5.m"

network_data_soc = PowerModels.parse_file(network_name)

PowerModels.print_summary(network_data_soc)

result_soc = run_opf(network_data_soc, SOCWRPowerModel, with_optimizer(Ipopt.Optimizer))

#result_soc = run_opf(network_data_soc, DCPPowerModel, with_optimizer(Ipopt.Optimizer))
#result_soc = run_opf(network_data_soc, QCWRPowerModel, with_optimizer(Ipopt.Optimizer))


network_data = PowerModels.parse_file(network_name)

result_soc

# Check for number of generators in the network
g_upper = length(network_data_soc["gen"])

for i=1:g_upper

    g=string(i)

    network_data["gen"][g]["pg"] = result_soc["solution"]["gen"][g]["pg"]

    network_data["gen"][g]["qg"] = result_soc["solution"]["gen"][g]["qg"]

end

result_ac = run_ac_pf(network_data, with_optimizer(Ipopt.Optimizer))

result_ac_unabr = run_ac_opf(network_data_soc, with_optimizer(Ipopt.Optimizer))

optimality_gap = 100*(result_ac_unabr["objective"] - result_soc["objective"])/result_ac_unabr["objective"] 
optimality_gap
result_ac

# Voltage magnitude accuracy
V_unabr = [bus["vm"] for (i,bus) in result_ac_unabr["solution"]["bus"]]
V_soc = sqrt.( [bus["w"] for (i,bus) in result_soc["solution"]["bus"]] )

v_err = V_unabr - V_soc
max_err = maximum( abs.(v_err) )
mean_err = sum( abs.(v_err) )/length(v_err)


 # Check for voltage violations

V_ac_sub = [bus["vm"] for (i,bus) in result_ac["solution"]["bus"]]

V_upper_lim = [bus["vmax"] for (i,bus) in network_data_soc["bus"]]

V_lower_lim = [bus["vmin"] for (i,bus) in network_data_soc["bus"]]


Check_lower = sum(V_ac_sub .< V_lower_lim)/ length(V_ac_sub)
Check_upper = sum(V_ac_sub .> V_upper_lim)/ length(V_ac_sub)

if Check_lower + Check_upper >0
    print("Constraint tightening needed")
else
    print("No constraint tightening needed")
end

# checking for flow violations







# Ranking most loaded lines







# Maximum phase angle deviation






# identifying active constraints which are getting violated





# Comparing solvers









 ## Constraint tightening

#  for i=1:length(network_data_soc["bus"])
#     g=string(i)

#     network_data_soc["bus"][g]["vmax"] = network_data_soc["bus"][g]["vmax"]*0.95
#     network_data_soc["bus"][g]["vmin"] = network_data_soc["bus"][g]["vmin"]*0.95

#  end

