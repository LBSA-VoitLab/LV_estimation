using Plots, Noise
using DifferentialEquations
using DelimitedFiles
plotlyjs()

noise_level = 0.05

t = 9
dts = .1

i1 = 1.
i2 = 1.
i3 = 1.
i4 = 1.

cc = -2
p1 = [2.5 cc]
p2 = [0.7 cc]
p3 = [0.2 cc]
p4 = [2.5 cc]
i12 = [-.5 .5]
i13 = [-1. 0.]
i14 = [-.5 -.5]
i23 = [1. 0.]
i24 = [.5 .5]
i34 = [0. 0.]
p12 = [p1 p2 i12]
p13 = [p1 p3 i13]
p14 = [p1 p4 i14]
p23 = [p2 p3 i23]
p24 = [p2 p4 i24]
p34 = [p3 p4 i34]
p123 = [p1 p2 p3 i12 i13 i23]
p124 = [p1 p2 p4 i12 i14 i24]
p134 = [p1 p3 p4 i13 i14 i34]
p234 = [p2 p3 p4 i23 i24 i34]
p1234 = [p1 p2 p3 p4 i12 i13 i14 i23 i24 i34]

tspan = (0.0,t)

function lv1(du,u,p,t)
  x1,x2 = u
  a1,b11 = p
  du[1] = dx1 = x1*(a1 + b11*x1)
  du[2] = dx2 = 0
end

u1 = [i1 0.]
u2 = [i2 0.]
u3 = [i3 0.]
u4 = [i4 0.]
u12 = [i1 i2]
u13 = [i1 i3]
u14 = [i1 i4]
u23 = [i1 i2]
u24 = [i2 i4]
u34 = [i3 i4]
u123 = [i1 i2 i3]
u124 = [i1 i2 i4]
u134 = [i1 i3 i4]
u234 = [i2 i3 i4]
u1234 = [i1 i2 i3 i4]

mon1 = ODEProblem(lv1,u1,tspan,p1)
mon2 = ODEProblem(lv1,u2,tspan,p2)
mon3 = ODEProblem(lv1,u3,tspan,p3)
mon4 = ODEProblem(lv1,u4,tspan,p4)

sol1 = DifferentialEquations.solve(mon1,SSPRK432(),saveat=1)
sol2 = DifferentialEquations.solve(mon2,SSPRK432(),saveat=1)
sol3 = DifferentialEquations.solve(mon3,SSPRK432(),saveat=1)
sol4 = DifferentialEquations.solve(mon4,SSPRK432(),saveat=1)
sol1 = add_gauss(sol1[:,2],noise_level)
sol2 = add_gauss(sol2[:,2],noise_level)
sol3 = add_gauss(sol3[:,2],noise_level)
sol4 = add_gauss(sol4[:,2],noise_level)
writedlm("m1.txt",[sol1.t sol1[1,:] sol1[2,:]])
writedlm("m2.txt",[sol2.t sol2[1,:] sol2[2,:]])
writedlm("m3.txt",[sol3.t sol3[1,:] sol3[2,:]])
writedlm("m4.txt",[sol4.t sol4[1,:] sol4[2,:]])


plot(sol1.t,sol1[1,:],grid=true,lw=3, label = "Bac1", title = "Growth Curves")
plot!(sol2.t,sol2[1,:], lw=3, label="Bac2")
plot!(sol3.t,sol3[1,:], lw=3, label="Bac3")
plot!(sol4.t,sol4[1,:], lw=3, label="Bac4")
current()

function lv2(du,u,p,t)
  x1,x2 = u
  a1,b11,a2,b22,b12,b21 = p
  du[1] = dx1 = x1*(a1 + b11*x1 + b12*x2)
  du[2] = dx2 = x2*(a2 + b21*x1 + b22*x2)
end

pw12 = ODEProblem(lv2,u12,tspan,p12)
pw13 = ODEProblem(lv2,u13,tspan,p13)
pw14 = ODEProblem(lv2,u14,tspan,p14)
pw23 = ODEProblem(lv2,u23,tspan,p23)
pw24 = ODEProblem(lv2,u24,tspan,p24)
pw34 = ODEProblem(lv2,u34,tspan,p34)

sol12 = DifferentialEquations.solve(pw12,SSPRK432(),saveat=1)
sol12 = add_gauss(sol12[:,2:3],noise_level)
plot(sol12.t,sol12[1,:],grid=true,lw=3, label = "Bac1", title = "Pairwise")
plot!(sol12.t,sol12[2,:], lw=3, label="Bac2")
writedlm("p12.txt",[sol12.t sol12[1,:] sol12[2,:]])

sol13 = DifferentialEquations.solve(pw13,SSPRK432(),saveat=1)
sol13 = add_gauss(sol13[:,2:3],noise_level)
plot(sol13.t,sol13[1,:],grid=true,lw=3, label = "Bac1", title = "Pairwise")
plot!(sol13.t,sol13[2,:], lw=3, label="Bac3")
writedlm("p13.txt",[sol13.t sol13[1,:] sol13[2,:]])

sol14 = DifferentialEquations.solve(pw14,SSPRK432(),saveat=1)
sol14 = add_gauss(sol14[:,2:3],noise_level)
plot(sol14.t,sol14[1,:],grid=true,lw=3, label = "Bac1", title = "Pairwise")
plot!(sol14.t,sol14[2,:], lw=3, label="Bac4")
writedlm("p14.txt",[sol14.t sol14[1,:] sol14[2,:]])

sol23 = DifferentialEquations.solve(pw23,SSPRK432(),saveat=1)
sol23 = add_gauss(sol23[:,2:3],noise_level)
plot(sol23.t,sol23[1,:],grid=true,lw=3, label = "Bac2", title = "Pairwise")
plot!(sol23.t,sol23[2,:], lw=3, label="Bac3")
writedlm("p23.txt",[sol23.t sol23[1,:] sol23[2,:]])

sol24 = DifferentialEquations.solve(pw24,SSPRK432(),saveat=1)
sol24 = add_gauss(sol24[:,2:3],noise_level)
plot(sol24.t,sol24[1,:],grid=true,lw=3, label = "Bac2", title = "Pairwise")
plot!(sol24.t,sol24[2,:], lw=3, label="Bac4")
writedlm("p24.txt",[sol24.t sol24[1,:] sol24[2,:]])

sol34 = DifferentialEquations.solve(pw34,SSPRK432(),saveat=1)
sol34 = add_gauss(sol34[:,2:3],noise_level)
plot(sol34.t,sol34[1,:],grid=true,lw=3, label = "Bac3", title = "Pairwise")
plot!(sol34.t,sol34[2,:], lw=3, label="Bac4")
writedlm("p34.txt",[sol34.t sol34[1,:] sol34[2,:]])

current()

function lv3(du,u,p,t)
  x1,x2,x3 = u
  a1,b11,a2,b22,a3,b33,b12,b21,b13,b31,b23,b32 = p
  du[1] = dx1 = x1*(a1 + b11*x1 + b12*x2 + b13*x3)
  du[2] = dx2 = x2*(a2 + b21*x1 + b22*x2 + b23*x3)
  du[3] = dx3 = x3*(a3 + b31*x1 + b32*x2 + b33*x3)
end

t123 = ODEProblem(lv3,u123,tspan,p123)
t124 = ODEProblem(lv3,u124,tspan,p124)
t134 = ODEProblem(lv3,u134,tspan,p134)
t234 = ODEProblem(lv3,u234,tspan,p234)

sol123 = DifferentialEquations.solve(t123,SSPRK432(),saveat=1)
sol123 = add_gauss(sol123[:,2:4],noise_level)
plot(sol123.t,sol123[1,:],grid=true,lw=3, label = "Bac1", title = "Triplicate")
plot!(sol123.t,sol123[2,:], lw=3, label="Bac2")
plot!(sol123.t,sol123[3,:], lw=3, label="Bac3")
writedlm("t123.txt",[sol123.t sol123[1,:] sol123[2,:] sol123[3,:]])

sol124 = DifferentialEquations.solve(t124,SSPRK432(),saveat=1)
sol124 = add_gauss(sol124[:,2:4],noise_level)
plot(sol124.t,sol124[1,:],grid=true,lw=3, label = "Bac1", title = "Triplicate")
plot!(sol124.t,sol124[2,:], lw=3, label="Bac2")
plot!(sol124.t,sol124[3,:], lw=3, label="Bac4")
writedlm("t124.txt",[sol124.t sol124[1,:] sol124[2,:] sol124[3,:]])

sol134 = DifferentialEquations.solve(t134,SSPRK432(),saveat=1)
sol134 = add_gauss(sol134[:,2:4],noise_level)
plot(sol134.t,sol134[1,:],grid=true,lw=3, label = "Bac1", title = "Triplicate")
plot!(sol134.t,sol134[2,:], lw=3, label="Bac3")
plot!(sol134.t,sol134[3,:], lw=3, label="Bac4")
writedlm("t134.txt",[sol134.t sol134[1,:] sol134[2,:] sol134[3,:]])

sol234 = DifferentialEquations.solve(t234,SSPRK432(),saveat=1)
sol234 = add_gauss(sol234[:,2:4],noise_level)
plot(sol234.t,sol234[1,:],grid=true,lw=3, label = "Bac2", title = "Triplicate")
plot!(sol234.t,sol234[2,:], lw=3, label="Bac3")
plot!(sol234.t,sol234[3,:], lw=3, label="Bac4")
writedlm("t234.txt",[sol234.t sol234[1,:] sol234[2,:] sol234[3,:]])

function lv4(du,u,p,t)
  x1,x2,x3,x4 = u
  a1,b11,a2,b22,a3,b33,a4,b44,b12,b21,b13,b31,b14,b41,b23,b32,b24,b42,b34,b43 = p
  du[1] = dx1 = x1*(a1 + b11*x1 + b12*x2 + b13*x3 + b14*x4)
  du[2] = dx2 = x2*(a2 + b21*x1 + b22*x2 + b23*x3 + b24*x4)
  du[3] = dx3 = x3*(a3 + b31*x1 + b32*x2 + b33*x3 + b34*x4)
  du[4] = dx4 = x4*(a4 + b41*x1 + b42*x2 + b43*x3 + b44*x4)
end

full = ODEProblem(lv4,u1234,tspan,p1234)

sol1234 = DifferentialEquations.solve(full,SSPRK432(),saveat=1)
sol1234 = add_gauss(sol1234[:,2:5],noise_level)
plot(sol1234.t,sol1234[1,:],grid=true,lw=3, label = "Bac1", title = "Full")
plot!(sol1234.t,sol1234[2,:], lw=3, label="Bac2")
plot!(sol1234.t,sol1234[3,:], lw=3, label="Bac3")
plot!(sol1234.t,sol1234[4,:], lw=3, label="Bac4")
writedlm("f.txt",[sol1234.t sol1234[1,:] sol1234[2,:] sol1234[3,:] sol1234[4,:]])
