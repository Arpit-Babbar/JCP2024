using Downloads: download
using Trixi, TrixiLW

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

initial_condition = initial_condition_constant

solver = DGSEM(polydeg=6, surface_flux=flux_lax_friedrichs,
               volume_integral = TrixiLW.VolumeIntegralFR(TrixiLW.LW()))


###############################################################################
# Get the uncurved mesh from a file (downloads the file if not available locally)

# Unstructured mesh with 48 cells of the square domain [-1, 1]^n
mesh_file = joinpath(@__DIR__, "GingerbreadManN6.inp")

# Map the unstructured mesh with the mapping above
mesh = P4estMesh{2}(mesh_file, initial_refinement_level=0)

semi = SemidiscretizationHyperbolic(mesh, get_time_discretization(solver),
                                    equations, initial_condition, solver,
                                    boundary_conditions=Dict(
                                      :Body => BoundaryConditionDirichlet(initial_condition),
                                      :Button1 => BoundaryConditionDirichlet(initial_condition),
                                      :Button2 => BoundaryConditionDirichlet(initial_condition),
                                      :Eye1 => BoundaryConditionDirichlet(initial_condition),
                                      :Smile => BoundaryConditionDirichlet(initial_condition),
                                      :Bowtie => BoundaryConditionDirichlet(initial_condition)
                                    ))


###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10.0)
lw_update = TrixiLW.semidiscretize(semi, get_time_discretization(solver), tspan);

summary_callback = SummaryCallback()

analysis_interval = 100000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=10000,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim,
				     output_directory="ginger")

callbacks = (summary_callback, analysis_callback, alive_callback, save_solution)

###############################################################################
# run the simulation

# run the simulation
time_int_tol = 1e-5 # Works best, it seems
tolerances = (; abstol=time_int_tol, reltol=time_int_tol);
dt_initial = 1e-6;
cfl_number = 0.3
sol = TrixiLW.solve_lwfr(lw_update, callbacks, dt_initial, tolerances,
#   time_step_computation = TrixiLW.Adaptive(),
  time_step_computation=TrixiLW.CFLBased(cfl_number),
#   limiters=(; stage_limiter!)
);
summary_callback() # print the timer summary
