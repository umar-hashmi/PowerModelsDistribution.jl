using Clp
using JuMP

struct Edge
    from; to; cost; capacity
end

edges =[Edge(1,2,1,0.5), Edge(1,3,2,0.4), Edge(1,4,3,0.6), Edge(2,5,1,0.5),Edge(4,5,2,0.5),Edge(3,5,2,0.6)]
mcf = Model(Clp.Optimizer)
@variable(mcf, 0 <= flow[e in edges] <= e.capacity)
@constraint(mcf, sum(flow[e] for e in edges if e.to==5) ==1)
@constraint(mcf, flowcon[node=2:4], sum(flow[e] for e in edges if e.to==node)  == sum(flow[e] for e in edges if e.from==node))
@objective(mcf, Min, sum(e.cost * flow[e] for e in edges))

print(mcf)

status =  optimize(mcf)