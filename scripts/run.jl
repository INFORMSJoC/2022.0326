@info "Loading packages..."
using MPI, Decomposition

Decomposition.setup_logger()

algorithm = ARGS[1]
instance = ARGS[2]
demand_scale = parse(Float64, ARGS[3])
limit_scale = parse(Float64, ARGS[4])
careful = parse(Bool, ARGS[5])

T = 24
max_time = 3600
max_iterations = 10000
verbose_solvers = [1]

if algorithm == "tcuc-central"
    @time @show solve_uc_centralized(instance,
                                     demand_scale=demand_scale,
                                     limit_scale=limit_scale,
                                     T=T,
                                     transmission=true,
                                     security=false)
elseif algorithm == "scuc-central"
    @time @show solve_uc_centralized(instance,
                                     demand_scale=demand_scale,
                                     limit_scale=limit_scale,
                                     transmission=true,
                                     security=true,
                                     T=T)
else
    MPI.Init()
    if algorithm == "scuc-isf"
        solve_uc_isf(instance,
                     demand_scale=demand_scale,
                     limit_scale=limit_scale,
                     T=T,
                     transmission=true,
                     security=true,
                     max_time=max_time,
                     max_iterations=max_iterations,
                     verbose_solvers=verbose_solvers,
                     careful=careful)
    elseif algorithm == "tcuc-isf"
        solve_uc_isf(instance,
                     demand_scale=demand_scale,
                     limit_scale=limit_scale,
                     T=T,
                     transmission=true,
                     security=false,
                     max_time=max_time,
                     max_iterations=max_iterations,
                     verbose_solvers=verbose_solvers)
    elseif algorithm == "tcuc-theta"
        solve_uc_theta(instance,
                       demand_scale=demand_scale,
                       limit_scale=limit_scale,
                       T=T,
                       max_time=max_time,
                       verbose_solvers=verbose_solvers)
   else
       @error("invalid algorithm: $algorithm")
   end
   MPI.Finalize()
end
