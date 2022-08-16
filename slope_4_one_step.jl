using DifferentialEquations, Plots, JuMP, Ipopt, DelimitedFiles
plotlyjs()

readin = readdlm("###FILE###.txt")


omit = 0

en = length(readin[:,1])

# Range
ran = 1:(en-omit)

t = readin[ran,1]
data1 = readin[ran,2]
data2 = readin[ran,3]
data3 = readin[ran,4]
data4 = readin[ran,5]

num = length(data1) - 1
b = zeros(num)
c = zeros(num)
d = zeros(num)
f = zeros(num)
a = zeros(num,5)

for i in 1:num
    a[i,1] = 1.0
    a[i,2] = data1[i]
    a[i,3] = data2[i]
    a[i,4] = data3[i]
    a[i,5] = data4[i]
    b[i] = (log(data1[i+1]) - log(data1[i]))/(t[i+1]-t[i])
    c[i] = (log(data2[i+1]) - log(data2[i]))/(t[i+1]-t[i])
    d[i] = (log(data3[i+1]) - log(data3[i]))/(t[i+1]-t[i])
    f[i] = (log(data4[i+1]) - log(data4[i]))/(t[i+1]-t[i])
end

reg = 0
p = zeros(4,5)

m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter"=>100000000,"print_level"=>2))
@variable(m,x[1:5])
@objective(m,Min,sum((a*x-b).^2) + reg*sum(x.^2))
@constraint(m,cap,x[2]<=0)
@constraint(m,growth,x[1]>=0)
@time status = JuMP.optimize!(m)
p1 = JuMP.value.(x)
er1 = sum((a*p1-b).^2)

m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter"=>100000,"print_level"=>2))
@variable(m,x[1:5])
@objective(m,Min,sum((a*x-c).^2) + reg*sum(x.^2))
@constraint(m,cap,x[3]<=0)
@constraint(m,growth,x[1]>=0)
@time status = JuMP.optimize!(m)
p2 = JuMP.value.(x)
er2 = sum((a*p2-c).^2)

m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter"=>100000,"print_level"=>2))
@variable(m,x[1:5])
@objective(m,Min,sum((a*x-d).^2) + reg*sum(x.^2))
@constraint(m,cap,x[4]<=0)
@constraint(m,growth,x[1]>=0)
@time status = JuMP.optimize!(m)
p3 = JuMP.value.(x)
er3 = sum((a*p3-d).^2)

m = Model(optimizer_with_attributes(Ipopt.Optimizer,"max_iter"=>100000,"print_level"=>2))
@variable(m,x[1:5])
@objective(m,Min,sum((a*x-f).^2) + reg*sum(x.^2))
@constraint(m,cap,x[5]<=0)
@constraint(m,growth,x[1]>=0)
@time status = JuMP.optimize!(m)
p4 = JuMP.value.(x)
er4 = sum((a*p4-f).^2)


p = [p1;p2;p3;p4]
println("")
println(p[1:5])
println(p[6:10])
println(p[11:15])
println(p[16:20])
println(er1+er2+er3+er4,"\t", "Error")
println([data1[end], data2[end], data3[end], data4[end]],"\t","Final Data Steady State")

function lotka_volterra(du,u,p,t)
  x1,x2,x3,x4 = u
  a1,b11,b12,b13,b14,a2,b21,b22,b23,b24,a3,b31,b32,b33,b34,a4,b41,b42,b43,b44 = p
  du[1] = dx1 = x1*(a1 + b11*x1 + b12*x2 + b13*x3 + b14*x4)
  du[2] = dx2 = x2*(a2 + b21*x1 + b22*x2 + b23*x3 + b24*x4)
  du[3] = dx3 = x3*(a3 + b31*x1 + b32*x2 + b33*x3 + b34*x4)
  du[4] = dx4 = x4*(a4 + b41*x1 + b42*x2 + b43*x3 + b44*x4)
end

u0 = [data1[1],data2[1],data3[1],data4[1]]
tspan = (0.0,t[end])
prob = ODEProblem(lotka_volterra,u0,tspan,p)
sol = DifferentialEquations.solve(prob,SSPRK432(),saveat=0.001,maxiters=1e7)
println(sol[end],"\t","Simulated Steady State")

plot(sol.t,sol[1,:],color=:blue, legend=:false, grid=:true)
plot!(sol.t,sol[2,:],color=:red)
plot!(sol.t,sol[3,:],color=:purple)
plot!(sol.t,sol[4,:],color=:orange)
scatter!(t,data1,color=:blue)
scatter!(t,data2,color=:red)
scatter!(t,data3,color=:purple)
scatter!(t,data4,color=:orange)
