using PowerModels
using Ipopt
using Plots
using JuMP

data = PowerModels.parse_file("data/matpower/case30.m")
ac_result = run_opf(data, ACPPowerModel, with_optimizer(Ipopt.Optimizer))
soc_result = run_opf(data, SOCWRPowerModel, with_optimizer(Ipopt.Optimizer))

# it is important to check the solver terminated nominally!
@assert  ac_result["termination_status"] == OPTIMAL ||  ac_result["termination_status"] == LOCALLY_SOLVED
@assert soc_result["termination_status"] == OPTIMAL || soc_result["termination_status"] == LOCALLY_SOLVED

# by convention optimality gaps are given as a percentage
gap = 100*(ac_result["objective"] - soc_result["objective"])/ac_result["objective"]

return (gap=gap, ub_time=ac_result["solve_time"], lb_time=soc_result["solve_time"]) 



