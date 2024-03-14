using PlotlyJS

begin # muss je nach Ordner angepasst werden
    include("C:/Users/liemann/Desktop/PowerDynamicsExtension/PowerDynamicsExtension.jl/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/Desktop/PowerDynamicsExtension/PowerDynamicsExtension.jl/src/grids/Example_Callbacks_OLTC.jl")
end

tspan  = (0.0,20.0)
pgsol = simulate_Example_OLTC(tspan);
plotallvoltages(pgsol);
plot([myplot(pgsol,"bus_load",:P0),myplot(pgsol,"bus_load",:Q0)])
