using Ipopt
using PowerModels
using CSV
using DelimitedFiles
using Distributions
using LinearAlgebra
using Plots

solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)
s = Dict("output" => Dict("branch_flows" => true) ) 

network_name = "data/matpower/case30.m"

network_data = PowerModels.parse_file(network_name)
network_solution = PowerModels.parse_file(network_name)

network_data_orig = PowerModels.parse_file(network_name)

result_orig = run_ac_opf(network_data_orig, with_optimizer(Ipopt.Optimizer); setting = s)



for (i,bus) in network_data["load"]
    network_data["load"][i]["pd"] = 1.7*deepcopy(network_data_orig["load"][i]["pd"] )
    network_data["load"][i]["qd"] = 1.7*deepcopy(network_data_orig["load"][i]["qd"])

    network_solution["load"][i]["pd"] = 1.7*deepcopy(network_data_orig["load"][i]["pd"] )
    network_solution["load"][i]["qd"] = 1.7*deepcopy(network_data_orig["load"][i]["qd"])
end

for (i,bus) in result_ac_original["solution"]["gen"]
    network_data["gen"][i]["pg"] = deepcopy(result_orig["solution"]["gen"][i]["pg"] )
    network_data["gen"][i]["qg"] = deepcopy(result_orig["solution"]["gen"][i]["qg"])

    network_solution["gen"][i]["pg"] = deepcopy(result_orig["solution"]["gen"][i]["pg"] )
    network_solution["gen"][i]["qg"] = deepcopy(result_orig["solution"]["gen"][i]["qg"])
end

for (i,bus) in result_ac_original["solution"]["bus"]
    network_data["bus"][i]["vm"] = deepcopy(result_orig["solution"]["bus"][i]["vm"] )
    network_solution["bus"][i]["vm"] = deepcopy(result_orig["solution"]["bus"][i]["vm"] )
end

#network_data["bus"]["2"]["bus_type"]= 3

network_solution["load"]["19"]["qd"] = network_solution["load"]["19"]["qd"] - 0.1
network_solution["load"]["19"]["pd"] = network_solution["load"]["19"]["pd"] - 0.2
network_solution["load"]["21"]["qd"] = network_solution["load"]["21"]["qd"] - 0.06
#network_solution["load"]["14"]["pd"] = network_solution["load"]["14"]["pd"] - 0.06

result_ac_inf = run_ac_pf(network_data, with_optimizer(Ipopt.Optimizer); setting = s)

result_solution = run_ac_pf(network_solution, with_optimizer(Ipopt.Optimizer); setting = s)



l_infla = PowerModels.component_table(network_data,"load", ["pd", "qd","load_bus"])
l_orig = PowerModels.component_table(network_data_orig,"load", ["pd", "qd"])
l_solution = PowerModels.component_table(network_solution,"load", ["pd", "qd"])

V_lim = PowerModels.component_table(network_data, "bus", ["vmin", "vmax"])
V_des = PowerModels.component_table(network_data, "bus", ["bus_type", "base_kv", "zone"])
Ang_lim = PowerModels.component_table(network_data, "branch", ["angmin", "angmax"])

V_infl = PowerModels.component_table(result_ac_inf["solution"], "bus", ["vm", "va"])
V_orig = PowerModels.component_table(result_orig["solution"], "bus", ["vm", "va"])
V_solution = PowerModels.component_table(result_solution["solution"], "bus", ["vm", "va"])

voltage_matrix = Array{Float64}([V_lim V_infl[:,2] V_orig[:,2] V_solution[:,2]])
voltage_angle = Array{Float64}([V_lim[:,1] V_infl[:,3] V_orig[:,3] V_solution[:,3]])
x_volt= voltage_matrix[:,1]
y_volt= voltage_matrix[:,2:end]
y_angle = voltage_angle[:,2:end]

plot(x_volt, y_volt, title = "Voltage Magnitude",label = ["Vmin" "Vmax" "AC infla" "AC OPF" "Solution"],legend=:bottomleft)
plot(x_volt, y_angle, title = "Voltage Angle",label = ["Min ang" "Max ang" "AC infla" "AC OPF" "Solution"],legend=:bottomright)

voltage_matrix_difference = Array{Float64}([V_lim[:,1] V_infl[:,2]-V_orig[:,2] V_solution[:,2]-V_orig[:,2]])
voltage_angle_difference = Array{Float64}([V_lim[:,1] V_infl[:,3]-V_orig[:,3] V_solution[:,3]-V_orig[:,3]])

plot(voltage_matrix_difference[:,1], voltage_matrix_difference[:,2:end], title = "Nominal voltage difference",label = ["AC infla" "V1" "V2" "V3" "V4" "V5" "V6"],legend=:bottomleft, legendfontsize=6)
plot(voltage_angle_difference[:,1], voltage_angle_difference[:,2:end], title = "Nominal voltage angle difference",label = ["AC infla" "V1" "V2" "V3" "V4" "V5" "V6"])


# bus type (1 = PQ, 2 = PV, 3 = ref, 4 = isolated)

############################################
# Branch analysis
############################################

PowerModels.update_data!(network_data, result_ac_inf["solution"])
dc_flows_v1 = PowerModels.calc_branch_flow_dc(network_data)
ac_flows_v1 = PowerModels.calc_branch_flow_ac(network_data)
ac_flows_limit = PowerModels.calc_branch_flow_ac(network_data)

PowerModels.update_data!(network_solution, result_solution["solution"])
ac_flows_solution = PowerModels.calc_branch_flow_ac(network_solution)

branc_details = PowerModels.component_table(network_data, "branch", ["f_bus", "t_bus", "rate_a","br_status","br_r", "br_x"])
from_bus_branch = branc_details[:,2]
to_bus_branch = branc_details[:,3]
R_branch = branc_details[:,6]
X_branch = branc_details[:,7]
Z_branch = R_branch + (X_branch)im
Z_branch_abs = abs.(Z_branch)*100/maximum(abs.(Z_branch))
Y_branch = transpose(1/Z_branch)
G_branch = real(Y_branch)  #Conductance
S_branch = imag(Y_branch)  #Suseptance

RbyXratio = sum(R_branch)/sum(X_branch)

# plot(R_branch./X_branch)

# plot(R_branch)
# plot!(X_branch)
# plot!(abs.(Z_branch))


Injection_ac_orig = PowerModels.component_table(result_orig["solution"],"branch",["qf", "pf", "qt","pt"])
Injection_ac_limit = PowerModels.component_table(ac_flows_limit,"branch",["qf", "pf", "qt","pt"])
Injection_ac_solution = PowerModels.component_table(ac_flows_solution,"branch",["qf", "pf", "qt","pt"])


thermal_ac_original = max( sqrt.(Injection_ac_orig[:,2].^2 + Injection_ac_orig[:,3].^2 ), sqrt.(Injection_ac_orig[:,4].^2 + Injection_ac_orig[:,5].^2 ))
thermal_ac_limit = max( sqrt.(Injection_ac_limit[:,2].^2 + Injection_ac_limit[:,3].^2 ), sqrt.(Injection_ac_limit[:,4].^2 + Injection_ac_limit[:,5].^2 ))
thermal_ac_solution = max( sqrt.(Injection_ac_solution[:,2].^2 + Injection_ac_solution[:,3].^2 ), sqrt.(Injection_ac_solution[:,4].^2 + Injection_ac_solution[:,5].^2 ))

losses_ac_orig = 100*abs.( Injection_ac_orig[:,3] + Injection_ac_orig[:,5] )./branc_details[:,4]
losses_ac_limit = 100*abs.( Injection_ac_limit[:,3] + Injection_ac_limit[:,5] )./branc_details[:,4]
losses_ac_solution = 100*abs.( Injection_ac_solution[:,3] + Injection_ac_solution[:,5] )./branc_details[:,4]

thr_loading_ac_orig = 100*thermal_ac_original./branc_details[:,4]
thr_loading_limit = 100*thermal_ac_limit./branc_details[:,4]
thr_loading_solution = 100*thermal_ac_solution./branc_details[:,4]

x_thermal = Injection_ac_orig[:,1]
y_thermal = [thr_loading_ac_orig thr_loading_limit thr_loading_solution]
y_loss = [losses_ac_orig losses_ac_limit losses_ac_solution]

plot(x_thermal, y_thermal, title = "thermal loading",label = ["AC OPF" "AC infla" "Solution" "V2" "V3" "V4" "V5" "V6"])

plot(x_thermal, y_loss, title = "line losses",label = ["AC OPF" "AC infla" "Solution" "V2" "V3" "V4" "V5" "V6"])




indx_max = max( maximum(from_bus_branch ), maximum(to_bus_branch ))
len_m = length(from_bus_branch)
branch_matrix = zeros(indx_max,indx_max)
rateA_matrix = zeros(indx_max,indx_max)
thermal_ac_matrix = zeros(indx_max,indx_max)
thermal_v1_matrix = zeros(indx_max,indx_max)
thermal_v6_matrix = zeros(indx_max,indx_max)

for id3 =1:len_m
    for id1=1:indx_max
        for id2 = 1: indx_max
            if id1 == branc_details[id3,2]  &&  id2 == branc_details[id3,3]
                branch_matrix[id1,id2] = Z_branch_abs[id3]
                branch_matrix[id2,id1] = Z_branch_abs[id3]

                rateA_matrix[id1,id2] = branc_details[id3,4]
                rateA_matrix[id2,id1] = branc_details[id3,4]

                thermal_ac_matrix[id1,id2] = thr_loading_ac_orig[id3]
                thermal_ac_matrix[id2,id1] = thr_loading_ac_orig[id3]

                thermal_v1_matrix[id1,id2] = thermal_ac_limit[id3]
                thermal_v1_matrix[id2,id1] = thermal_ac_limit[id3]
            end
        end
    end
end

pyplot()
# plt = plot(branch_matrix, st=:heatmap, clim=(0,100), color=:coolwarm, colorbar_title="Impedances of branch", xlabel = "to bus", ylabel = "from bus")
# #savefig(plt, "branchImpedance.eps")
# plt2 = plot(rateA_matrix, st=:heatmap, color=:coolwarm, colorbar_title="RateA of network branches", xlabel = "to bus", ylabel = "from bus")
# #savefig(plt2, "rateAnetwork.eps")
# plt3 = plot(thermal_ac_matrix, st=:heatmap, color=:coolwarm, colorbar_title="AC thermal", xlabel = "to bus", ylabel = "from bus")
# #savefig(plt3, "thermal_nominal.eps")
# plt4 = plot(thermal_v1_matrix, st=:heatmap, color=:coolwarm, colorbar_title="V1 thermal", xlabel = "to bus", ylabel = "from bus")
# #savefig(plt4, "thermal_v1.eps")

##################################################
# Generator analysis
##################################################

Generator_limits = PowerModels.component_table(network_data,"gen",["gen_bus", "pmax", "pmin", "qmax","qmin"])
Generator_output = PowerModels.component_table(result_orig["solution"],"gen",["pg", "qg"])
Generator_output_infla = PowerModels.component_table(result_ac_inf["solution"],"gen",["pg", "qg"])
Generator_output_solution = PowerModels.component_table(result_solution["solution"],"gen",["pg", "qg"])

x_gen = Generator_limits[:,1]
y_gen_act = [Generator_limits[:,3] Generator_limits[:,4] Generator_output[:,2] Generator_output_infla[:,2] Generator_output_solution[:,2] ]
#y_gen_act_percent = 

y_gen_REact = [Generator_limits[:,5] Generator_limits[:,6] Generator_output[:,3] Generator_output_infla[:,3] Generator_output_solution[:,3]]

plot(x_gen, y_gen_act, title = "Generator active power",label = ["Pmax" "Pmin" "AC OPF" "AC infla" "solution" "V2" "V3" "V4" "V5" "V6"], lw=[3 3 1 1 1 1 1 1 1 1])

plot(x_gen, y_gen_REact, title = "Generator reactive power",label = ["Qmax" "Qmin" "AC OPF" "AC infla" "solution" "V2" "V3" "V4" "V5" "V6"], lw=[3 3 1 1 1 1 1 1 1 1], legendfontsize=6,legend=:bottomright)
