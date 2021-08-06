using GLMakie; GLMakie.activate!()
using DifferentialEquations

###################################################################################
### Model from https://journals.aps.org/pre/abstract/10.1103/PhysRevE.99.052205 ###
###################################################################################


const Cm = 0.5 #uF/cm^2 membrane capacitance
const gCa = 0.025 #mS/cm^2 max conductance Ca2+ channels
const gK = 0.04 #mS/cm^2 max conductance K+ channels
const VCa = 100.0 #mV Ca2+ reversal potential
const VK = -80.0 #mV K+ reversal potential
const taun = 300.0 #ms K+ activation time constant
const tauh = 80.0 #ms Ca+ activation time constant
const Vm = -35.0 #mV voltage at midpoint minf
const sm = 6.24 #mV slope param for minf
const Vn = -40.0 #mV voltage at midpoint ninf
const sn = 5.0 #mV slope param for ninf
const Vh = -20.0 #mV voltage at midpoint sinf
const sh = 8.6 #mV slope parameter sinf

function minf(V)
    1/(1+exp((Vm-V)/sm))
end

function ninf(V)
    1/(1+exp((Vn-V)/sn))
end

function hinf(V)
    1/(1+exp((V-Vh)/sh))
end

function Isti(t,PCL)
    k = floor(Int, t/PCL)
    40.0*Float64(k*PCL < t < k*PCL + 1.0)
end

function vo_luo_rudy!(du, u, p, t)
    PCL = p
    V, n, h = u
    IK = gK*n*(V-VK)
    ICa = gCa*minf(V)*h*(V-VCa)
    Istiv = Isti(t, PCL)
    du[1] = dV = (-(IK + ICa) + Istiv)/Cm
    du[2] = dn = (ninf(V) - n)/taun
    du[3] = dh = (hinf(V) - h)/tauh
end

###############################################################

function tr(PCL)
    u0 = zeros(3)
    tspan = (0.0,20000.0)
    p = PCL

    prob = ODEProblem(vo_luo_rudy!, u0, tspan, p)
    sol = solve(prob,RK4(),adaptive=false,dt=.2);

    sol
end

fig = Figure()
ax = Axis(fig[1,1])
sl = Slider(fig[2, 1], range = 1100:.05:1600, startvalue = 1100.0)

obs = sl.value
sol = @lift(tr($obs))

empty!(ax)

on(sol) do sol
    empty!(ax)
    lines!(ax,sol.t, sol[1,:],color=:blue)
    #println(obs.val)
end
