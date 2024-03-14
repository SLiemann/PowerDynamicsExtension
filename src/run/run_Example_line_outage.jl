using PlotlyJS

begin # muss je nach Ordner angepasst werden
    include("C:/Users/liemann/Desktop/PowerDynamicsExtension/PowerDynamicsExtension.jl/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/Desktop/PowerDynamicsExtension/PowerDynamicsExtension.jl/src/grids/Example_Line_Outage.jl")
end

tspan  = (0.0,8.0)
pgsol = simulate_line_outage(tspan);
plotallvoltages(pgsol);
plot(myplot(pgsol,"bus_sm",:θ))
plot(myplot(pgsol,"bus_sm",:ω))

plot([myplot(pgsol,"bus_load",:P0),myplot(pgsol,"bus_load",:Q0)])
