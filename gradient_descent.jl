using DifferentialEquations, DiffEqSensitivity, Plots, Flux, DiffEqFlux, DelimitedFiles, OptimizationFlux, OptimizationOptimisers, Optimization

# File location
file = "###.txt"

# Approximate carrying capacity of the system (e.g. 1e10)
cc = 1e10

# What is your initial parameter guess? Use zeros(# of params to estimate) if unsure
p0 = zeros(20)

# Number of iterations
its = 1e6

# Penalty term
λ = 0.0

# No penalty
penalty(x) = 0

# LASSO
# penalty(x) = λ*sum(abs,x)

# Ridge
# penalty(x) = λ*x^2

# Elastic Net
# penalty(x) = λ*sum(abs,x) + λ*x^2




### Begin code
a = readdlm(file,Float64)
omit = 0

t = a[1:end-omit,1]'
data1 = a[1:end-omit,2]'*scale
data2 = a[1:end-omit,3]'*scale
data3 = a[1:end-omit,4]'*scale
data4 = a[1:end-omit,5]'*scale
data0 = [data1;data2;data3;data4]

function lotka_volterra!(du,u,p,t)
  x1,x2,x3,x4 = u
  a1,b11,b12,b13,b14,a2,b21,b22,b23,b24,a3,b31,b32,b33,b34,a4,b41,b42,b43,b44 = p
  du[1] = dx1 = x1*(a1 + b11*x1 + b12*x2 + b13*x3 + b14*x4)
  du[2] = dx2 = x2*(a2 + b21*x1 + b22*x2 + b23*x3 + b24*x4)
  du[3] = dx3 = x3*(a3 + b31*x1 + b32*x2 + b33*x3 + b34*x4)
  du[4] = dx4 = x4*(a4 + b41*x1 + b42*x2 + b43*x3 + b44*x4)
end


u0 = [data1[1],data2[1],data3[1],data4[1]]
t_span = (t[1], t[end])
prob0 = DiffEqBase.ODEProblem(lotka_volterra!, u0, t_span, p0)

passages = []
condition(u,t,integrator) = t ∈ passages
affect!(integrator) = integrator.u /= 40
cb = DiscreteCallback(condition,affect!)

function loss_function(p)
  try
    tmp_prob = remake(prob0,p=p)
    prediction = solve(tmp_prob, Tsit5(), callback=cb, tstops = passages, saveat = t)
    loss = sum(abs2,Array(prediction) - data0) + sum(penalty,p)
    return loss
  catch
    return 1e50
  end
end

optf = OptimizationFunction((x,p) -> loss_function(x), Optimization.AutoForwardDiff())
prob = Optimization.OptimizationProblem(optf, p0)
res = solve(prob, OptimizationOptimisers.AdaMax(1e-3), maxiters = its, progress = true)

p1 = res.u
p1[2] = p1[2]*scale
p1[3] = p1[3]*scale
p1[4] = p1[4]*scale
p1[5] = p1[5]*scale
p1[7] = p1[7]*scale
p1[8] = p1[8]*scale
p1[9] = p1[9]*scale
p1[10] = p1[10]*scale
p1[12] = p1[12]*scale
p1[13] = p1[13]*scale
p1[14] = p1[14]*scale
p1[15] = p1[15]*scale
p1[17] = p1[17]*scale
p1[18] = p1[18]*scale
p1[19] = p1[19]*scale
p1[20] = p1[20]*scale
u0 = u0/scale

prob1 = DiffEqBase.ODEProblem(lotka_volterra!, u0, t_span, res.u)
sol1 = solve(prob1, Tsit5(), callback=cb, tstops = passages, saveat = 0.01)
plot(sol1.t,sol1[1,:],color=[1], legend=:true,label="x1 fit", grid=:true, lw=3, yaxis=:log10)
plot!(sol1.t,sol1[2,:],color=[2], label="x2 fit", lw=3)
plot!(sol1.t,sol1[3,:],color=[3], label="x3 fit", lw=3)
plot!(sol1.t,sol1[4,:],color=[4], label="x4 fit", lw=3)
scatter!(t, data1/scale, color = [1], legend =:false, label = "bac1")
scatter!(t, data2/scale, color = [2], label = "bac2")
scatter!(t, data3/scale, color = [3], label = "bac3")
scatter!(t, data4/scale, color = [4], label = "bac4")