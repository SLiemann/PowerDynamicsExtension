using PlotlyJS, DataFrames
using ModelingToolkit   

begin
    include("C:/Users/liemann/Desktop/PowerDynamicsExtension/PowerDynamicsExtension.jl/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/Desktop/PowerDynamicsExtension/PowerDynamicsExtension.jl/src/grids/Example_Callbacks_OLTC.jl")
end

tspan  = (0.0,25.0)
pgsol = simulate_Example_OLTC(tspan);
plotallvoltages(pgsol);
