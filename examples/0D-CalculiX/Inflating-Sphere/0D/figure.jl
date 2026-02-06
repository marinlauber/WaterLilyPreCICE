using CSV,DataFrames,Plots
# load the file
# data = CSV.read(joinpath(@__DIR__,"ConstantFlow.csv"), DataFrame; delim = ',')
# data = CSV.read(joinpath(@__DIR__,"Windkessel_open.csv"), DataFrame; delim = ',')
data = CSV.read(joinpath(@__DIR__,"Windkessel_closed.csv"), DataFrame; delim = ',')
# find all the index of the last iteration of each time step
idx_last = [maximum([findall(data.timestep.≈i)...,1]) for i in 1:maximum(Int,data.timestep)]
# make new dataframe with only the last iteration of each time step
data_conv = @views(data[idx_last, :])

# plot the number of iterations per timestep
lsim = maximum(data.timestep)

# maximum number of iterations
xp = collect(1:maximum(data.iter))

# all the first iteration values index
time = data_conv.time
deltat = time[2]-time[1]
tmax = maximum(time)

p1 = plot(time, data_conv.Plv_0D, label="Pv 0D", lw=3, c=:cyan, ylabel="Pressure (mmHg)", alpha=0.5,
           xlims=(0,tmax), ylims=(0,80))
plot!(p1, time, data_conv.Pact_3D, label="Pact 3D", lw=3, ls=:dot, c=3, alpha=0.5)
plot!(p1, time, data_conv.Plv_3D, label="Plv 3D", lw=3, c=1)
plot!(p1, time, data_conv.Plv_0D, label="Plv 0D", lw=3, ls=:dashdotdot, c=2)
plot!(p1, time, data_conv.Pao_0D, label="Pa 0D", lw=3, ls=:dashdotdot, c=:blue)

p6 = plot(data_conv.Vlv_3D, data_conv.Plv_3D, label="3D", ylabel="Plv (mmHg)",
          lw=3, xlabel="Vlv (ml)", xlims=(0,160), ylims=(0,80))
plot!(p6, data_conv.Vlv_0D, data_conv.Plv_0D, label="0D", ls=:dashdotdot, lw=3, c=:orange)

# make the combined plot
plot(p1, p6, layout=(1,2), size=(1200,600), dpi=300)
savefig("pressure_volume_sphere.png")