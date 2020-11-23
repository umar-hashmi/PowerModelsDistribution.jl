# Test script to use the infiniteOPF package

# Optimisation problem: series connected R-L-C ceircuit
# u(t) = ur(t) + ul(t) + uc(t)
# obj:  min ∫u(t)dt
# s.t.  u(t) = r*i(t) + l*di(t)/dt + c*∫i(t)dt
#       0.2 ≤ i(t) ≤ 1A
#       0.5 ≤ u(t) ≤ 1V

# Declare packages: 
# You will have to add them, e.g. Pkg.add("Ipopt") etc.
using InfiniteOpt, JuMP, Distributions, Ipopt

# instantiate model
m = InfiniteModel(Ipopt.Optimizer)

# define parameters
T = 10  # time range 0-10s
r = 1 #1 Ohm
l = 1e-3 # 1 mH
c = 1e-6 # 1 μF

# define time parameter
@infinite_parameter(m, t ∈ [0, T], num_supports = 10)

# define variables
@infinite_variable(m, 0.5 ≤ u(t) ≤ 1)
@infinite_variable(m, 0.2 ≤ i(t) ≤ 1)

# define objective
@objective(m, Min, integral(u, t))

# define circuit constraint
@constraint(m, r*i + l * @deriv(i, t) + c * @integral(i, t) == u)

#solve optimisation
optimize!(m)

#print result
opt_objective = objective_value(m)