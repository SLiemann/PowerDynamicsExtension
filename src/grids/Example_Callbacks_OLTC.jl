using PowerDynamics
using OrderedCollections: OrderedDict
import DiffEqBase: initialize_dae!

function pg_cb_oltc() 
    buses=OrderedDict(
        "bus0" => SlackAlgebraic(U=1.0), 
        "bus1" => VoltageDependentLoad(P=0.0, Q=0.0, U=1.0, A=1.0, B=0.0,Y_n = complex(0.0)),
        "bus_load" => GeneralVoltageDependentLoad(P=-0.75, Q = -0.3, U=1.0, Ap=0.0, Bp=1.0,Aq = 1.0, Bq= 0.00,Y_n = complex(0.0)))

    branches=OrderedDict(
        "Line_0-1"=> PiModelLine(from= "bus0", to = "bus1",y=1.0/(0.05+1im*0.1), y_shunt_km=0.0, y_shunt_mk=0.0),
        "OLTC"=> StaticPowerTransformerTapParam(from="bus1",to="bus_load",Sbase=100e6,Srated=100e6,uk=0.10,XR_ratio=Inf,
                                           i0=0.0,Pv0=0.0,tap_side = "LV",tap_pos = 0,tap_inc = 1.0,tap_max=20,tap_min=-20,p_ind=1))
        return PowerGrid(buses, branches)
end

function Initialize_pg_cb_oltc()
    pg = pg_cb_oltc()
    U1,δ1,ic0,cu = PowerFlowClassic(pg,iwamoto = false,max_tol = 1e-4)
    pg, ic = InitializeInternalDynamics(pg,ic0)
    display("V => δ")
    display("------")
    display(U1.=>δ1)
    return pg,ic
end

function simulate_Example_OLTC(tspan::Tuple{Float64,Float64})
    pg,ic  = Initialize_pg_cb_oltc()
    params = [0] # tap position of OLT 
    problem = ODEProblem{true}(rhs(pg),ic,tspan,params)
   
    timer_start = -1.0
    tap_dir = 1.0

    branch_oltc = "OLTC"
    index_U_oltc = PowerDynamics.variable_index(pg.nodes,pg.lines[branch_oltc].to,:u_r) # at LV side, otherwise "from"

    function TapState(integrator)
        timer_start = integrator.t
        integrator.p[1] += 1*tap_dir #position of parameter have to be selected automatically 
        initialize_dae!(integrator,BrownFullBasicInit())
        auto_dt_reset!(integrator)
    end

    function voltage_deadband(u,t,integrator)
         0.99 <= sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) <= 1.01 # deadband and other information could be stored in node 
    end

    function timer_off(integrator)
        if timer_start != -1
            timer_start = -1
        end
    end

    function voltage_outside_low(u,t,integrator)
         sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) < 0.99
    end

    function voltage_outside_high(u,t,integrator)
        sqrt(u[index_U_oltc]*u[index_U_oltc] + u[index_U_oltc+1]*u[index_U_oltc+1]) > 1.01
   end

    function timer_on_low(integrator)
        tap_dir = 1
        if timer_start == -1
            timer_start = integrator.t
        end
    end

    function timer_on_high(integrator)
        tap_dir = -1
        if timer_start == -1
            timer_start = integrator.t
        end
    end

    function timer_hit(u,t,integrator)
        if timer_start == -1
            return false
        else
            return t-timer_start >= 1.0 # waiting time of OLTC until tap change when voltage is outside deadband 
        end
    end

    cb1 = DiscreteCallback(voltage_deadband, timer_off) 
    cb2 = DiscreteCallback(voltage_outside_low, timer_on_low)
    cb3 = DiscreteCallback(voltage_outside_high, timer_on_high)
    cb4 = DiscreteCallback(timer_hit, TapState)
    
    sol = solve(problem, Rodas4(autodiff=true), callback = CallbackSet(cb1,cb2,cb3,cb4), dtmax = 1e-3,force_dtmin=false,maxiters=1e6, initializealg = BrownFullBasicInit(),alg_hints=:stiff,abstol=1e-9,reltol=1e-9) 
    return PowerGridSolution(sol, pg)
end