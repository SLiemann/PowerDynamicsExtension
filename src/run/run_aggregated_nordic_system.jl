using PlotlyJS, DataFrames
using ModelingToolkit   

Sbase = 8000e6
Ubase = 400e3
Zbase = Ubase^2/Sbase

begin
    include("C:/Users/liemann/Desktop/PowerDynamicsExtension/PowerDynamicsExtension.jl/src/include_custom_nodes_lines_utilities.jl")
    include("C:/Users/liemann/Desktop/PowerDynamicsExtension/PowerDynamicsExtension.jl/src/grids/Aggregated_Nordic_Test_System.jl")
    nothing
end

pg0,ic0 = Initialize_N32_GFM_TS();
@time pgsol = simulate_LTVS_N32_simulation_TS(pg0,ic0,(0.0,30.0),(20.0+1im*20)/Zbase);
p1 = plotallvoltages(pgsol);
plot(myplot(pgsol,"bus_gfm",:i_abs))
plot(myplot(pgsol,"bus_gfm",:q_imax))
plot(myplot(pgsol,"bus_gfm",:q_idcmax))
plot(myplot(pgsol,"bus_gfm",:idc0))
plot([myplot(pgsol,"bus_gfm",:θ),myplot(pgsol1,"bus_gfm",:θ)])
plot([myplot(pgsol,"bus_gfm",:P0),myplot(pgsol0_per,"bus_gfm",:P0)])

