using PowerDynamics
using OrderedCollections: OrderedDict
import DiffEqBase: initialize_dae!

function pg_line_outage() 
    buses=OrderedDict(
        "bus0" => SlackAlgebraic(U=1.0), 
        "bus_fault" => ThreePhaseFault(rfault=8e3,xfault=8e3,p_ind=collect(1:2)),
        "bus1" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_sm" => gentpjAVROEL(Sbase=100e6,Srated=100e6, H=6.0, P=0.8, D=0.0, Ω=50, R_a=0, T_d0s=7.0, T_q0s=1.5, T_d0ss=0.05, T_q0ss=0.05, X_d=2.2, X_q=2.0, X_ds=0.3, X_qs=0.4, X_dss=0.2, X_qss=0.2, X_l=0.15, S_10=0.1, S_12=0.3,K_is=0.0, V0 = 1.0, Ifdlim = 3.0618, L1 = -20.0, G1 = 120.0, Ta = 5.0, Tb = 12.5, G2 = 10.0, L2 = 5.0),
        "bus_load" => GeneralVoltageDependentLoad(P=-0.5, Q = -0.2, U=1.0, Ap=0.0, Bp=1.0,Aq = 1.0, Bq= 0.00,Y_n = complex(0.0)))

    branches=OrderedDict(
        "Line_0-f"=> PiModelLine(from= "bus0", to = "bus_fault",y=1.0/(0.005+1im*0.1), y_shunt_km=0.0, y_shunt_mk=0.0),
        "Line_f_sm"=> PiModelLine(from= "bus_fault", to = "bus_sm",y=1.0/(0.005+1im*0.1), y_shunt_km=0.0, y_shunt_mk=0.0),
        "Line_0-1"=> PiModelLine(from= "bus0", to = "bus1",y=1.0/(0.005+1im*0.1), y_shunt_km=0.0, y_shunt_mk=0.0),
        "Line_1-sm"=> PiModelLine(from= "bus1", to = "bus_sm",y=1.0/(0.005+1im*0.2), y_shunt_km=0.0, y_shunt_mk=0.0),
        "Trafo"=> StaticPowerTransformerTapParam(from="bus1",to="bus_load",Sbase=100e6,Srated=100e6,uk=0.10,XR_ratio=Inf, i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = 0,tap_inc = 1.0,tap_max=20,tap_min=-20,p_ind=3))
        return PowerGrid(buses, branches)
end

function Initialize_pg_line_outage()
    pg = pg_line_outage()
    Qmax   = [Inf,Inf,sqrt(1-0.8^2),Inf] #maxium reactive power of generator
    Qmin   = -Qmax
    U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-4,iter_max = 100,Qmax = Qmax, Qmin = Qmin,Qlimit_iter_check=5)
    pg, ic = InitializeInternalDynamics(pg,ic0)
    display("V => δ")
    display("------")
    display(U1.=>δ1)
    return pg,ic
end

function get_pg_after_line_outage(pg::PowerGrid)
    nodes_postfault = deepcopy(pg.nodes)
    branches_postfault = deepcopy(pg.lines)
    delete!(branches_postfault,"Line_f_sm")
    return PowerGrid(nodes_postfault,branches_postfault)
end

function get_continous_fault(pg::PowerGrid)
    pg_fault = deepcopy(pg)
    pg_fault.nodes["bus_fault"] = ThreePhaseFaultContinouos(rfault=8e3,xfault=8e3,Tf=1e-3/2,p_ind=collect(1:2))
    return pg_fault
end

function simulate_line_outage(tspan::Tuple{Float64,Float64})
    pg,ic  = Initialize_pg_line_outage()
    params = [1e4,1e4,0] # fault impedance and tap position of OLTC 
    problem = ODEProblem{true}(rhs(pg),ic,tspan,params)

    rfault = 0.005
    xfault = 0.005
    pg_postfault = get_pg_after_line_outage(pg)
   
    function short_circuit(integrator)        
        integrator.p[1] = rfault #rfault
        integrator.p[2] = xfault #xfault

        pg_cfault = get_continous_fault(pg);
        ic_init= deepcopy(integrator.sol[end])
        len = length(symbolsof(pg.nodes["bus1"]))
        len += length(symbolsof(pg.nodes["bus_sm"]))
        len += length(symbolsof(pg.nodes["bus_load"]))
        # insert two extra states for continouos fault
        ic_tmp = vcat(ic_init[1:end-len],[pg_cfault.nodes["bus_fault"].rfault,pg_cfault.nodes["bus_fault"].xfault],ic_init[end-len+1:end])
        # create problem and simulate for 10ms
        op_prob = ODEProblem(rhs(pg_cfault), ic_tmp, (0.0, 0.001), integrator.p)
        sol_tmp = solve(op_prob,Rodas4(),dtmax=1e-4,initializealg = BrownFullBasicInit(),alg_hints=:stiff,verbose=false,abstol=1e-8,reltol=1e-8)

        ic_end = sol_tmp.u[end]
        # delete states of continouos fault
        ic_end = vcat(ic_end[1:end-len-2],ic_end[end-len+1:end])
        # change only algebraic states of original problem
        ind_as = findall(x-> iszero(x),diag(integrator.f.mass_matrix))
        #ind_as = getVoltageSymbolPositions(pg)
        for i in ind_as
            ic_init[i] = ic_end[i]
        end
        
        integrator.u = deepcopy(ic_init)
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
    end

    function line_outage(integrator)
        integrator.p[1] = 1e4 #fault is zero again
        integrator.p[2] = 1e4 #fault is zero again

        #First create continouos fault and then post-fault grid
        pg_pcfault = get_continous_fault(pg);
        #pg_pcfault = GetPostFaultLTVSPG_TS(pg_pcfault);
        ic_init= deepcopy(integrator.sol[end])
        len = length(symbolsof(pg.nodes["bus1"]))
        len += length(symbolsof(pg.nodes["bus_sm"]))
        len += length(symbolsof(pg.nodes["bus_load"]))
        # insert two extra states for continouos fault
        ic_tmp = vcat(ic_init[1:end-len],[rfault,xfault],ic_init[end-len+1:end])
        # create problem and simulate for 20 ms
        op_prob = ODEProblem(rhs(pg_pcfault), ic_tmp, (0.0, 0.001), integrator.p)
        sol_tmp = solve(op_prob,Rodas4(),dtmax=1e-4,initializealg = BrownFullBasicInit(),alg_hints=:stiff,verbose=false,abstol=1e-8,reltol=1e-8)

        ic_end = sol_tmp.u[end]
        # delete states of continouos fault
        ic_end = vcat(ic_end[1:end-len-2],ic_end[end-len+1:end])
        # change only algebraic states of original problem
        ind_as = findall(x-> iszero(x),diag(integrator.f.mass_matrix))
        for i in ind_as
            ic_init[i] = ic_end[i]
        end

        integrator.u = deepcopy(ic_init)
        ode   = rhs(pg_postfault)
        integrator.f = ode
        integrator.cache.tf.f = integrator.f
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
    end

    cb1 = PresetTimeCallback(1.0, short_circuit)
    cb2 = PresetTimeCallback(1.15, line_outage)
    
    sol = solve(problem, Rodas4(autodiff=true), callback = CallbackSet(cb1,cb2), dtmax = 1e-3,force_dtmin=false,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-9,reltol=1e-9) 
    return PowerGridSolution(sol, pg)
end