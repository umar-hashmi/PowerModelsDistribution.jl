using Ipopt
using PowerModels
using SortingLab
using CSV
using DelimitedFiles
#using PyPlot
using Distributions

using LinearAlgebra
using Plots



solver = optimizer_with_attributes(Ipopt.Optimizer, "tol"=>1e-6, "print_level"=>0)

network_name = "data/matpower/case30.m"
#network_name = "C:\\Users\\uhashmi\\.julia\\packages\\PowerModels\\4FzNJ\\test\\data\\matpower\\case30.m"
#network_name = "data/matpower/case33bw.m"
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

for (i,bus) in result_ac_original["solution"]["bus"]
    network_data_new["bus"][i]["vm"] = deepcopy(result_ac_original["solution"]["bus"][i]["vm"] )
    network_data_new_v1["bus"][i]["vm"] = deepcopy(result_ac_original["solution"]["bus"][i]["vm"] )
    network_data_new_v2["bus"][i]["vm"] = deepcopy(result_ac_original["solution"]["bus"][i]["vm"] )
    network_data_new_v3["bus"][i]["vm"] = deepcopy(result_ac_original["solution"]["bus"][i]["vm"] )
    network_data_new_v4["bus"][i]["vm"] = deepcopy(result_ac_original["solution"]["bus"][i]["vm"] )
    network_data_new_v5["bus"][i]["vm"] = deepcopy(result_ac_original["solution"]["bus"][i]["vm"] )
    network_data_new_v6["bus"][i]["vm"] = deepcopy(result_ac_original["solution"]["bus"][i]["vm"] )
end

a=[1.1 1.2 1.1 1.05 0.95 0.9 0.8]
r = [1 1 1 1 1 1 1 1]
#r=[0.9 1.1 1.05 1.02 0.98 0.95 0.9]

for (i,bus) in network_data["load"]
    network_data_new["load"][i]["pd"] = a[1]*deepcopy( network_data["load"][i]["pd"] )
    network_data_new["load"][i]["qd"] = r[1]*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v1["load"][i]["pd"] = a[2]*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v1["load"][i]["qd"] = r[2]*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v2["load"][i]["pd"] = a[3]*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v2["load"][i]["qd"] = r[3]*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v3["load"][i]["pd"] = a[4]*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v3["load"][i]["qd"] = r[4]*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v4["load"][i]["pd"] = a[5]*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v4["load"][i]["qd"] = r[5]*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v5["load"][i]["pd"] = a[6]*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v5["load"][i]["qd"] = r[6]*deepcopy(network_data["load"][i]["qd"])

    network_data_new_v6["load"][i]["pd"] = a[7]*deepcopy( network_data["load"][i]["pd"] )
    network_data_new_v6["load"][i]["qd"] = r[7]*deepcopy(network_data["load"][i]["qd"])
end

# Analyzing load

l0 = PowerModels.component_table(network_data,"load", ["pd", "qd","load_bus"])
l1 = PowerModels.component_table(network_data_new,"load", ["pd", "qd"])
l2 = PowerModels.component_table(network_data_new_v1,"load", ["pd", "qd"])
l3 = PowerModels.component_table(network_data_new_v2,"load", ["pd", "qd"])
l4 = PowerModels.component_table(network_data_new_v3,"load", ["pd", "qd"])
l5 = PowerModels.component_table(network_data_new_v4,"load", ["pd", "qd"])
l6 = PowerModels.component_table(network_data_new_v5,"load", ["pd", "qd"])
l7 = PowerModels.component_table(network_data_new_v6,"load", ["pd", "qd"])

load_active_matrix = Array{Float64}([l0[:,4] l0[:,2] l1[:,2] l2[:,2] l3[:,2] l4[:,2] l5[:,2] l6[:,2] l7[:,2]])
load_REactive_matrix = Array{Float64}([l0[:,4] l0[:,3] l1[:,3] l2[:,3] l3[:,3] l4[:,3] l5[:,3] l6[:,3] l7[:,3]])
x_load= load_active_matrix[:,1]
y_active= load_active_matrix[:,2:end]
y_REactive= load_REactive_matrix[:,2:end]

plot(x_load, y_active , seriestype = :scatter, title = "Active power",label = ["AC OPF" "AC infla" "V1" "V2" "V3" "V4" "V5" "V6"])
plot(x_load, y_REactive , seriestype = :scatter, title = "Reactive power",label = ["AC OPF" "AC infla" "V1" "V2" "V3" "V4" "V5" "V6"])

active_matrix_difference = Array{Float64}([l0[:,4] l0[:,2]-l1[:,2] l0[:,2]-l2[:,2] l0[:,2]-l3[:,2] l0[:,2]-l4[:,2] l0[:,2]-l5[:,2] l0[:,2]-l6[:,2] l0[:,2]-l7[:,2]])
reactive_difference = Array{Float64}([l0[:,4] l0[:,3]-l1[:,3] l0[:,3]-l2[:,3] l0[:,3]-l3[:,3] l0[:,3]-l4[:,3] l0[:,3]-l5[:,3] l0[:,3]-l6[:,3] l0[:,3]-l7[:,3]])

plot(active_matrix_difference[:,1], active_matrix_difference[:,2:end], title = "Active power difference from nominal",label = ["AC infla" "V1" "V2" "V3" "V4" "V5" "V6"])

plot(active_matrix_difference[:,1], reactive_difference[:,2:end], title = "Reactive power difference from nominal",label = ["AC infla" "V1" "V2" "V3" "V4" "V5" "V6"])

# Power flow solutions

result_new = run_ac_opf(network_data_new, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v1 = run_ac_pf(network_data_new_v1, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v2 = run_ac_pf(network_data_new_v2, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v3 = run_ac_pf(network_data_new_v3, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v4 = run_ac_pf(network_data_new_v4, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v5 = run_ac_pf(network_data_new_v5, with_optimizer(Ipopt.Optimizer); setting = s)
result_new_v6 = run_ac_pf(network_data_new_v6, with_optimizer(Ipopt.Optimizer); setting = s)





####################################
# Bus analysis
####################################

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
y_angle = voltage_angle[:,2:end]
plot(x_volt, y_volt, title = "Voltage Magnitude",label = ["Vmin" "Vmax" "AC OPF" "AC infla" "V1" "V2" "V3" "V4" "V5" "V6"],legend=:bottomleft)
plot(x_volt, y_angle, title = "Voltage Angle",label = ["Min ang" "Max ang" "AC OPF" "AC infla" "V1" "V2" "V3" "V4" "V5" "V6"],legend=:bottomright)

voltage_matrix_difference = Array{Float64}([m1[:,1] m0[:,2]-m01[:,2] m0[:,2]-m2[:,2] m0[:,2]-m3[:,2] m0[:,2]-m4[:,2] m0[:,2]-m5[:,2] m0[:,2]-m6[:,2] m0[:,2]-m7[:,2]])
voltage_angle_difference = Array{Float64}([m1[:,1] m0[:,3]-m01[:,3] m0[:,3]-m2[:,3] m0[:,3]-m3[:,3] m0[:,3]-m4[:,3] m0[:,3]-m5[:,3] m0[:,3]-m6[:,3] m0[:,3]-m7[:,3]])

plot(voltage_matrix_difference[:,1], voltage_matrix_difference[:,2:end], title = "Nominal voltage difference",label = ["AC infla" "V1" "V2" "V3" "V4" "V5" "V6"],legend=:bottomleft, legendfontsize=6)

plot(voltage_angle_difference[:,1], voltage_angle_difference[:,2:end], title = "Nominal voltage angle difference",label = ["AC infla" "V1" "V2" "V3" "V4" "V5" "V6"])
xlabel!("Bus ID")

V_bar = 0.5.*(m1[:,2] + m1[:,3]);

del_voltage_matrix = Array{Float64}([m1[:,1] abs.(100*(m0[:,2]-V_bar)./V_bar) abs.(100*(m01[:,2]-V_bar)./V_bar) abs.(100*(m2[:,2]-V_bar)./V_bar) abs.(100*(m3[:,2]-V_bar)./V_bar) abs.(100*(m4[:,2]-V_bar)./V_bar) abs.(100*(m5[:,2]-V_bar)./V_bar) abs.(100*(m6[:,2]-V_bar)./V_bar) abs.(100*(m7[:,2]-V_bar)./V_bar)])





############################################
# Branch analysis
############################################

PowerModels.update_data!(network_data_new_v1, result_new_v1["solution"])
PowerModels.update_data!(network_data_new_v2, result_new_v2["solution"])
PowerModels.update_data!(network_data_new_v3, result_new_v3["solution"])
PowerModels.update_data!(network_data_new_v4, result_new_v4["solution"])
PowerModels.update_data!(network_data_new_v5, result_new_v5["solution"])
PowerModels.update_data!(network_data_new_v6, result_new_v6["solution"])

dc_flows_v1 = PowerModels.calc_branch_flow_dc(network_data_new_v1)
ac_flows_v1 = PowerModels.calc_branch_flow_ac(network_data_new_v1)
ac_flows_v2 = PowerModels.calc_branch_flow_ac(network_data_new_v2)
ac_flows_v3 = PowerModels.calc_branch_flow_ac(network_data_new_v3)
ac_flows_v4 = PowerModels.calc_branch_flow_ac(network_data_new_v4)
ac_flows_v5 = PowerModels.calc_branch_flow_ac(network_data_new_v5)
ac_flows_v6 = PowerModels.calc_branch_flow_ac(network_data_new_v6)
ac_flows_limit = PowerModels.calc_branch_flow_ac(network_data_new)

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

plot(R_branch./X_branch)

plot(R_branch)
plot!(X_branch)
plot!(abs.(Z_branch))


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

plot(x_thermal, y_thermal, title = "thermal loading",label = ["AC OPF" "AC infla" "V1" "V2" "V3" "V4" "V5" "V6"])

plot(x_thermal, y_loss, title = "line losses",label = ["AC OPF" "AC infla" "V1" "V2" "V3" "V4" "V5" "V6"])




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

                thermal_v1_matrix[id1,id2] = thr_loading_v1[id3]
                thermal_v1_matrix[id2,id1] = thr_loading_v1[id3]

                thermal_v6_matrix[id1,id2] = thr_loading_v6[id3]
                thermal_v6_matrix[id2,id1] = thr_loading_v6[id3]
            end
        end
    end
end

pyplot()
plt = plot(branch_matrix, st=:heatmap, clim=(0,100), color=:coolwarm, colorbar_title="Impedances of branch", xlabel = "to bus", ylabel = "from bus")
savefig(plt, "branchImpedance.eps")
plt2 = plot(rateA_matrix, st=:heatmap, color=:coolwarm, colorbar_title="RateA of network branches", xlabel = "to bus", ylabel = "from bus")
savefig(plt2, "rateAnetwork.eps")
plt3 = plot(thermal_ac_matrix, st=:heatmap, color=:coolwarm, colorbar_title="AC thermal", xlabel = "to bus", ylabel = "from bus")
savefig(plt3, "thermal_nominal.eps")
plt4 = plot(thermal_v1_matrix, st=:heatmap, clim=(0,100), color=:coolwarm, colorbar_title="V1 thermal", xlabel = "to bus", ylabel = "from bus")
savefig(plt4, "thermal_v1.eps")
plt5 = plot(thermal_v6_matrix, st=:heatmap, clim=(0,100), color=:coolwarm, colorbar_title="V6 thermal", xlabel = "to bus", ylabel = "from bus")
savefig(plt5, "thermal_v6.eps")



##################################################
# Generator analysis
##################################################

Generator_limits = PowerModels.component_table(network_data,"gen",["gen_bus", "pmax", "pmin", "qmax","qmin"])
Generator_output = PowerModels.component_table(result_ac_original["solution"],"gen",["pg", "qg"])
Generator_output_infla = PowerModels.component_table(result_new["solution"],"gen",["pg", "qg"])
Generator_output_v1= PowerModels.component_table(result_new_v1["solution"],"gen",["pg", "qg"])
Generator_output_v2= PowerModels.component_table(result_new_v2["solution"],"gen",["pg", "qg"])
Generator_output_v3= PowerModels.component_table(result_new_v3["solution"],"gen",["pg", "qg"])
Generator_output_v4= PowerModels.component_table(result_new_v4["solution"],"gen",["pg", "qg"])
Generator_output_v5= PowerModels.component_table(result_new_v5["solution"],"gen",["pg", "qg"])
Generator_output_v6= PowerModels.component_table(result_new_v6["solution"],"gen",["pg", "qg"])


x_gen = Generator_limits[:,1]
y_gen_act = [Generator_limits[:,3] Generator_limits[:,4] Generator_output[:,2] Generator_output_infla[:,2] Generator_output_v1[:,2] Generator_output_v2[:,2] Generator_output_v3[:,2] Generator_output_v4[:,2] Generator_output_v5[:,2] Generator_output_v6[:,2]]
#y_gen_act_percent = 

y_gen_REact = [Generator_limits[:,5] Generator_limits[:,6] Generator_output[:,3] Generator_output_infla[:,3] Generator_output_v1[:,3] Generator_output_v2[:,3] Generator_output_v3[:,3] Generator_output_v4[:,3] Generator_output_v5[:,3] Generator_output_v6[:,3]]

plot(x_gen, y_gen_act, title = "Generator active power",label = ["Pmax" "Pmin" "AC OPF" "AC infla" "V1" "V2" "V3" "V4" "V5" "V6"], lw=[3 3 1 1 1 1 1 1 1 1])

plot(x_gen, y_gen_REact, title = "Generator reactive power",label = ["Qmax" "Qmin" "AC OPF" "AC infla" "V1" "V2" "V3" "V4" "V5" "V6"], lw=[3 3 1 1 1 1 1 1 1 1], legendfontsize=6,legend=:bottomright)



#######################################
#### Voltage sensitivity matrix with active power #######
#######################################

using StatsPlots

l_bus = length(x_volt)

Voltage_sensitivity_matrix = zeros(l_bus,l_bus)
V_slack_plus = zeros(l_bus,l_bus)
V_slack_minus = zeros(l_bus,l_bus)
V_slack = zeros(l_bus,l_bus)

Nom_voltage_mat = PowerModels.component_table(result_ac_original["solution"], "bus", ["vm", "va"] )

Voltage_reference = (Nom_voltage_mat[:,2])

Voltage_lims = PowerModels.component_table(network_data, "bus", ["vmin", "vmax", "vm", "va"])

l0 = PowerModels.component_table(network_data,"load", ["pd", "qd","load_bus"])

#ld= PowerModels.component_table(network_duplicate,"load", ["pd", "qd","load_bus"])

bus_id_load = l0[:,4]
id_bus_row1 = l0[:,1]

indx_max = max( maximum(from_bus_branch ), maximum(to_bus_branch ))
len_m = length(from_bus_branch)

thermal_duplicate = zeros(indx_max,indx_max)
losses_duplicate = zeros(indx_max,indx_max)
thermal_slack = zeros(indx_max,indx_max)

thermal_accumulate = zeros(indx_max,indx_max)
losses_accumulate = zeros(indx_max,indx_max)

del_power = 0.1

for indx_s = 1:length(bus_id_load)

    rst = string(id_bus_row1[indx_s])
    index_v = bus_id_load[indx_s]

    network_duplicate = PowerModels.parse_file(network_name)

    for (i,bus) in result_ac_original["solution"]["gen"]
        network_duplicate["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
        network_duplicate["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])
    end

   network_duplicate["load"][rst]["pd"] = network_duplicate["load"][rst]["pd"] - del_power
 #  network_duplicate["load"][rst]["qd"] = network_duplicate["load"][rst]["qd"] - del_power

    #- min(0.1,network_duplicate["load"][rst]["pd"] )

    result_duplicate = run_ac_pf(network_duplicate, with_optimizer(Ipopt.Optimizer); setting = s)

    Duplicate_voltage_mat = PowerModels.component_table(result_duplicate["solution"], "bus", ["vm", "va"] )

    Voltage_dublicate = (Duplicate_voltage_mat[:,2])

    intermediate_matrix = transpose(Voltage_dublicate) - transpose(Voltage_reference)

    V_slack_plus[index_v,:] = transpose(Voltage_lims[:,3]) - transpose(Voltage_dublicate)
    V_slack_minus[index_v,:] = transpose(Voltage_dublicate) - transpose(Voltage_lims[:,2])
    V_slack[index_v,:] = min.(V_slack_plus[index_v,:], V_slack_minus[index_v,:])

    Voltage_sensitivity_matrix[index_v,:] = intermediate_matrix/del_power

    # branc_details = PowerModels.component_table(network_data, "branch", ["f_bus", "t_bus", "rate_a","br_status","br_r", "br_x"])

    ac_flows_duplicate = PowerModels.calc_branch_flow_ac(network_duplicate )
    Injection_ac_duplicate = PowerModels.component_table(ac_flows_duplicate,"branch",["qf", "pf", "qt","pt"])
    thermal_ac_duplicate = max( sqrt.(Injection_ac_duplicate[:,2].^2 + Injection_ac_duplicate[:,3].^2 ), sqrt.(Injection_ac_duplicate[:,4].^2 + Injection_ac_duplicate[:,5].^2 ))
    losses_ac_duplicate = 100*abs.( Injection_ac_duplicate[:,3] + Injection_ac_duplicate[:,5] )./branc_details[:,4]
    thr_loading_duplicate = 100*thermal_ac_duplicate./branc_details[:,4]

    for id3 =1:len_m
        for id1=1:indx_max
            for id2 = 1: indx_max
                if id1 == branc_details[id3,2]  &&  id2 == branc_details[id3,3]
                    thermal_duplicate[id1,id2] = thr_loading_duplicate[id3]
                    thermal_duplicate[id2,id1] = thr_loading_duplicate[id3]
                    thermal_slack[id1,id2] = 100 - thr_loading_duplicate[id3]
                    thermal_slack[id2,id1] = 100 - thr_loading_duplicate[id3]

                    losses_duplicate[id1,id2] = losses_ac_duplicate[id3]
                    losses_duplicate[id2,id1] = losses_ac_duplicate[id3]
                end
            end
        end
    end

    thermal_accumulate = thermal_accumulate + thermal_duplicate
    losses_accumulate = losses_accumulate + losses_duplicate

    # println(Voltage_sensitivity_matrix )
    # println(index_v)
    # println(indx_s)

end

heatmap(Voltage_sensitivity_matrix)

interpret_this = sum(Voltage_sensitivity_matrix,dims = 2)
interpret_this2 = sum(abs.(Voltage_sensitivity_matrix),dims = 2)

groupedbar(sum(Voltage_sensitivity_matrix,dims = 2), bar_position = :stack, bar_width=0.7, label="Cumulative effect", xlabel="bus id", ylabel="sum of voltage fluctions")

plt_rank = groupedbar(sum(abs.(Voltage_sensitivity_matrix),dims = 2), bar_position = :stack, bar_width=0.7, label="Cumulative effect", xlabel="bus id", ylabel="sum of |voltage fluctions|" )
savefig(plt_rank, "rank2.eps")

mutual_coupling_voltage = sum(abs.(Voltage_sensitivity_matrix) - abs.(Diagonal(Voltage_sensitivity_matrix)),dims = 2)

plt_mutual=groupedbar(mutual_coupling_voltage, bar_position = :stack, bar_width=0.7, label="Cumulative mutual effect", xlabel="bus id", ylabel="sum of |voltage fluctions|")
savefig(plt_mutual, "mutual2.eps")
plt_box = boxplot(transpose(Voltage_sensitivity_matrix), leg = false, xlabel="bus id", ylabel="voltage fluctuation", title="Voltage sensitivity")
savefig(plt_box, "vsens2.eps")

plt_slack_v = boxplot(transpose(V_slack), leg = false, xlabel="bus id", ylabel="V_slack", title="Voltage slackness")
savefig(plt_slack_v, "slackV2.eps")

plt_slack_thermal = boxplot(transpose(thermal_slack), leg = false, xlabel="bus id", ylabel="T_slack", title="Thermal slackness")
#savefig(plt_slack_thermal, "slackT.eps")

#heatmap(losses_accumulate/length(bus_id_load))

#heatmap(thermal_accumulate/length(bus_id_load))


plot(thermal_accumulate/length(bus_id_load), st=:heatmap,clim=(0,100), color=:coolwarm, colorbar_title="thermal aggregate loading", xlabel="bus id to", ylabel="bus id from")
plot(losses_accumulate/length(bus_id_load), st=:heatmap,clim=(0,5), color=:coolwarm, colorbar_title="losses aggregate", xlabel="bus id to", ylabel="bus id from")

plot(Diagonal(Voltage_sensitivity_matrix), st=:heatmap, color=:coolwarm, colorbar_title="losses aggregate")

plot(V_slack, st=:heatmap, color=:coolwarm, colorbar_title="V_slack")
plot(thermal_slack, st=:heatmap, color=:coolwarm, colorbar_title="thermal_slack", xlabel="bus id to", ylabel="bus id from")

plot(Voltage_sensitivity_matrix, st=:heatmap, color=:coolwarm, colorbar_title="Voltage sensitivity")

plt_self = heatmap(Diagonal(Voltage_sensitivity_matrix), xlabel="bus id", ylabel="bus id", title= "Self coupling",colorbar_title="Voltage sensitivity")
savefig(plt_self, "selfActive.eps")




writedlm( "FileName2.csv",  Voltage_sensitivity_matrix, ',')






#######################################
#### Voltage sensitivity matrix with reactive power #######
#######################################



Voltage_sensitivity_matrix_rea = zeros(30,30)

bus_id_load = l0[:,4]
id_bus_row1 = l0[:,1]

for indx_s = 1:length(bus_id_load)

    rst = string(id_bus_row1[indx_s])
    index_v = bus_id_load[indx_s]

    network_duplicate = PowerModels.parse_file(network_name);

    for (i,bus) in result_ac_original["solution"]["gen"]
        network_duplicate["gen"][i]["pg"] = deepcopy(result_ac_original["solution"]["gen"][i]["pg"] )
        network_duplicate["gen"][i]["qg"] = deepcopy(result_ac_original["solution"]["gen"][i]["qg"])
    end

    network_duplicate["load"][rst]["qd"] = network_duplicate["load"][rst]["qd"] + 0.1

    result_duplicate = run_ac_pf(network_duplicate, with_optimizer(Ipopt.Optimizer); setting = s)

    Duplicate_voltage_mat = PowerModels.component_table(result_duplicate["solution"], "bus", ["vm", "va"] )

    Voltage_dublicate = (Duplicate_voltage_mat[:,2])

    intermediate_matrix = transpose(Voltage_dublicate) - transpose(Voltage_reference)

    Voltage_sensitivity_matrix_rea[index_v,:] = intermediate_matrix

    # println(Voltage_sensitivity_matrix )
    # println(index_v)
    # println(indx_s)

end

heatmap(Voltage_sensitivity_matrix_rea)

groupedbar(sum(Voltage_sensitivity_matrix_rea,dims = 2), bar_position = :stack, bar_width=0.7)




### Load profile

load_pf = l0[:,2]./sqrt.(l0[:,2].^2 + l0[:,3].^2)

plot(l0[:,4], load_pf)
plot!(l0[:,4], l0[:,3])
plot!(l0[:,4], l0[:,2])

mean_pf = sum(l0[:,2])/sqrt(sum(l0[:,2])^2 + sum(l0[:,3])^2)
P_tot = sum(l0[:,2])
Q_tot = sum(l0[:,3])
ϵ = 0.3
ζ = 0.4
