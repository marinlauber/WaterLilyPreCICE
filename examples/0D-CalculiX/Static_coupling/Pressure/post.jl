using CSV,DataFrames,Plots

let
    step_data = CSV.read("/home/marin/Workspace/WaterLilyPreCICE/examples/0D-CalculiX/Static_coupling/Pressure/sphere_output.csv", DataFrame; delim = ',')
    #start plots
    p1 = plot()
    p2 = plot()
    p4 = plot()
    p5 = plot()

    # plot the number of iterations per timestep
    lsim = maximum(step_data.timestep)
    iterations = diff([searchsortedfirst(step_data.timestep, i) for i in 1:lsim])
    p3 = plot(1:lsim-1, iterations, label=:none, xlabel="time step", ylabel="iterations",
              ylims=(0,ceil(1.2maximum(iterations))))

    # maximum number of iterations
    xp = collect(1:maximum(step_data.iter))

    # loop over results
    for i in 1:maximum(step_data.timestep)
        # find all iteration during this time step
        idx_ts = findall(step_data.timestep.≈i)
        # plot between t and t+Δt
        xs = range(i-1,i,length=length(idx_ts)+1)[1:end-1]
        # volume and pressure change raw
        plot!(p1, xs, step_data.pressure[idx_ts]*100, label=:none, c=:red)
        plot!(p2, xs, step_data.volume[idx_ts], label=:none, c=:black)
        # convergence of the volume
        plot!(p4, xs, abs.(step_data.target_volume[idx_ts] .- step_data.volume[idx_ts])./step_data.target_volume[idx_ts],
              label=:none, c=:black, yscale=:log10, xlabel="time step", ylabel="relative volume residual",
              xlims=(0,10),ylims=(1e-6,1))
        # flow rate
        last_v = maximum([findall(step_data.timestep.≈i-1)...,1])
        plot!(p5, xs, step_data.volume[idx_ts] .- step_data.volume[last_v], label=:none, c=:blue,
              xlabel="time step", ylabel="Flow rate",xlims=(0,10), ylims=(0,.1))
    end
    plot!(p1, [0], [0], c=:red, label=:none, ylabel="Pressure", xlims=(0,10), ylims=(0,1.5))
    plot!(p2, [0], [0], c=:black, label=:none, ylabel="Volume", xlims=(0,10), ylims=(0,1.5))

    # make the combined plot
    plot(p1, p2, p3, p4, p5, layout=(2,3), size=(1200,600))
    savefig("results.png")
end
