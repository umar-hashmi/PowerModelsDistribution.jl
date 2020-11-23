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

#network_name = "data/matpower/pglib_opf_case39_epri.m"
network_name = "data/matpower/pglib_opf_case57_ieee.m"


s = Dict("output" => Dict("branch_flows" => true) ) 

network_data_soc = PowerModels.parse_file(network_name)

PowerModels.print_summary(network_data_soc)

#result_soc = run_opf(network_data_soc, SOCWRPowerModel, with_optimizer(Ipopt.Optimizer); setting = s)

result_soc = run_opf(network_data_soc, DCPPowerModel, with_optimizer(Ipopt.Optimizer); setting = s)
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

result_ac = run_ac_pf(network_data, with_optimizer(Ipopt.Optimizer); setting = s)

result_ac_unabr = run_ac_opf(network_data_soc, with_optimizer(Ipopt.Optimizer); setting = s)

optimality_gap = 100*(result_ac_unabr["objective"] - result_soc["objective"])/result_ac_unabr["objective"] 
optimality_gap
result_ac

# Voltage magnitude accuracy
V_unabr = [bus["vm"] for (i,bus) in result_ac_unabr["solution"]["bus"]]

V_soc = [bus["vm"] for (i,bus) in result_soc["solution"]["bus"]]    # use this for DC approx
#V_soc = sqrt.( [bus["w"] for (i,bus) in result_soc["solution"]["bus"]] )       # use this for soc 

v_err = V_unabr - V_soc
max_err = maximum( abs.(v_err) )
mean_err = sum( abs.(v_err) )/length(v_err)


# Check for voltage violations
V_ac_sub = [bus["vm"] for (i,bus) in result_ac["solution"]["bus"]]
l_sub = length(V_ac_sub)
V_upper_lim = [bus["vmax"] for (i,bus) in network_data_soc["bus"]]
l_up = length(V_upper_lim)
V_lower_lim = [bus["vmin"] for (i,bus) in network_data_soc["bus"]]
l_low = length(V_lower_lim)
v_length = length(result_ac["solution"]["bus"])

#for j=1:v_length
#    g=string(j)
    
if l_sub == l_up && l_sub == l_low

    Check_lower = sum(V_ac_sub .< V_lower_lim)/ length(V_ac_sub)
    Check_upper = sum(V_ac_sub .> V_upper_lim)/ length(V_ac_sub)

    Check_lower_soc = sum(V_soc .< V_lower_lim)/ length(V_soc)
    Check_upper_soc = sum(V_soc .> V_upper_lim)/ length(V_soc)


    Check_lower_unabr = sum(V_unabr .< V_lower_lim)/ length(V_unabr)
    Check_upper_unabr = sum(V_unabr .> V_upper_lim)/ length(V_unabr)

else
    V_lower_lim = V_lower_lim[1:l_sub]
    V_upper_lim = V_upper_lim[1:l_sub]

    Check_lower = sum(V_ac_sub .< V_lower_lim)/ length(V_ac_sub)
    Check_upper = sum(V_ac_sub .> V_upper_lim)/ length(V_ac_sub)

    Check_lower_soc = sum(V_soc .< V_lower_lim)/ length(V_soc)
    Check_upper_soc = sum(V_soc .> V_upper_lim)/ length(V_soc)


    Check_lower_unabr = sum(V_unabr .< V_lower_lim)/ length(V_unabr)
    Check_upper_unabr = sum(V_unabr .> V_upper_lim)/ length(V_unabr)
end


if Check_lower + Check_upper >0
    print("Constraint tightening needed")
else
    print("No constraint tightening needed")
end

if Check_lower_unabr + Check_upper_unabr >0
    print("Constraint tightening needed in unabridged case")
else
    print("No constraint tightening needed in unabridged case")
end

# Individual voltage violations

plot(V_lower_lim)
plot!(V_upper_lim)
plot!(V_ac_sub)
plot!(V_soc)
plot!(V_unabr)

result_table = [result_ac_unabr["objective"] Check_lower_unabr+Check_upper_unabr result_soc["objective"] optimality_gap mean_err max_err Check_lower+Check_upper]

# checking for flow violations

thermal_limit= sort( [branch["rate_a"] for (i,branch) in network_data_soc["branch"]] )

thermal_soc = sort(abs.( [branch["pt"] for (i,branch) in result_soc["solution"]["branch"]])  )

thermal_ac_unab  = sort( abs.( [branch["pt"] for (i,branch) in result_ac_unabr["solution"]["branch"]] ) )

loading_soc = 100*thermal_soc./thermal_limit

loading_ac_unab = 100*thermal_ac_unab./thermal_limit


r11 = result_ac_unabr["objective"]
r12 = Check_lower_unabr+Check_upper_unabr
r13 = maximum(loading_ac_unab)

r21 = result_soc["objective"]
r22 = optimality_gap
r23 = mean_err
r24 = max_err
r25 = Check_lower+Check_upper
r26 = maximum(loading_soc)

result_table = [ r11  r12  r13  r21  r22  r23  r24  r25  r26]


plot(loading_ac_unab)
plot!(loading_soc)
#plot!(thermal_limit)



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

