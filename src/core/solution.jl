function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:JuMP.VariableRef,1})
    return JuMP.value.(var.data)
end


function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:JuMP.NonlinearExpression,1})
    return JuMP.value.(var.data)
end


function _IM.build_solution_values(var::JuMP.Containers.DenseAxisArray{<:Number,1})
    return var.data
end
