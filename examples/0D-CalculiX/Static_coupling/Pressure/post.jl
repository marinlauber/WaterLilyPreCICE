using CSV,DataFrames,Plots
# load the file
step_data = CSV.read("/home/marin/Workspace/WaterLilyPreCICE/examples/0D-CalculiX/Static_coupling/Pressure/sphere_output.csv", DataFrame; delim = ',')
let
    #start plots
    p4 = plot()
    p5 = plot()

    # plot the number of iterations per timestep
    lsim = maximum(step_data.timestep)
    iterations = diff([searchsortedfirst(step_data.timestep, i) for i in 1:lsim])
    p3 = plot(1:lsim-1, iterations, label=:none, xlabel="time step", ylabel="iterations",
              ylims=(0,ceil(1.2maximum(iterations))))

    # maximum number of iterations
    xp = collect(1:maximum(step_data.iter))

    # all the first iteration values index
    idx_1 = [maximum([findall(step_data.timestep.≈i)...,1]) for i in 1:maximum(Int,step_data.timestep)]
    time = step_data.time[idx_1]
    deltat = time[2]-time[1]
    tmax = maximum(time)
    p6 = plot(time, step_data.VLV[idx_1], label="VLV", alpha=0.5, c=:orange, xlims=(0,tmax))
    plot!(p6, time, step_data.PAO[idx_1], label="Pao", alpha=0.5, c=:blue)
    
    p1 = plot(time, step_data.pressure[idx_1], c=:red, label=:none, ylabel="Pressure (mmHg)", alpha=0.5, xlims=(0,tmax))
    p2 = plot(time, step_data.volume[idx_1], c=:black, label=:none, ylabel="Volume (ml)", alpha=0.5, xlims=(0,tmax), ylims=(0,160))
    # plot!(p1, time, step_data.Plv[idx_1]./6000, label="Plv", alpha=0.5, c=:green)

    # flow rates
    p5 = plot(time, -step_data.Qao[idx_1]./12, label=:"Qao",
              xlabel="time (s)", ylabel="Flow rate (ml/s)",xlims=(0,tmax))
    plot!(p5, time, step_data.Qmv[idx_1]./12, label="Qmv")

    # loop over results
    for i in 1:maximum(Int,step_data.timestep)
        # find all iteration during this time step
        idx_ts = findall(step_data.timestep.≈i)
        # plot between t and t+Δt
        xs = range(time[i]-deltat,time[i],length=length(idx_ts)+1)[2:end]
        # volume and pressure change raw
        plot!(p1, xs, step_data.pressure[idx_ts], label=:none, c=:red, lw=2)
        plot!(p2, xs, step_data.volume[idx_ts], label=:none, c=:black, lw=2)
        # convergence of the volume
        plot!(p4, xs, abs.(step_data.target_volume[idx_ts] .- step_data.volume[idx_ts])./step_data.target_volume[idx_ts],
              label=:none, c=:black, yscale=:log10, xlabel="time (s)", ylabel="relative volume residual",
              xlims=(0,tmax),ylims=(1e-6,1))
        # flow rate using the previous time step last volume (converged)
        last_v = maximum([findall(step_data.timestep.≈i-1)...,1])
        plot!(p5, xs, step_data.volume[idx_ts] .- step_data.volume[last_v], label=:none, c=:blue, lw=2)
        # 0D model ouptu
        plot!(p6, xs, step_data.VLV[idx_ts], label=:none, c=:orange, lw=2)
        plot!(p6, xs, step_data.PAO[idx_ts], label=:none, c=:blue, lw=2)
        # plot!(p1, xs, step_data.Plv[idx_ts]./60000, label=:none, c=:green, lw=2)
    end

    # make the combined plot
    plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1200,600))
    # savefig("results.png")
end
