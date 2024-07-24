using Downloads: download
using Trixi, TrixiLW

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

initial_condition = initial_condition_constant

solver = DGSEM(polydeg=5, surface_flux=flux_lax_friedrichs,
               volume_integral = TrixiLW.VolumeIntegralFR(TrixiLW.LW()))


###############################################################################
# Get the uncurved mesh from a file (downloads the file if not available locally)

# Mapping as described in https://arxiv.org/abs/2012.12040 but reduced to 2D
function mapping(xi_, eta_)
  # Transform input variables between -1 and 1 onto [0,3]
  xi = 1.5 * xi_ + 1.5
  eta = 1.5 * eta_ + 1.5

  y = eta + 3/8 * (cos(1.5 * pi * (2 * xi - 3)/3) *
                   cos(0.5 * pi * (2 * eta - 3)/3))

  x = xi + 3/8 * (cos(0.5 * pi * (2 * xi - 3)/3) *
                  cos(2 * pi * (2 * y - 3)/3))

  return SVector(x, y)
end

###############################################################################
# Get the uncurved mesh from a file (downloads the file if not available locally)

# Unstructured mesh with 48 cells of the square domain [-1, 1]^n
mesh_file = joinpath(@__DIR__, "square_unstructured_1.inp")
isfile(mesh_file) || download("https://gist.githubusercontent.com/efaulhaber/a075f8ec39a67fa9fad8f6f84342cbca/raw/a7206a02ed3a5d3cadacd8d9694ac154f9151db7/square_unstructured_1.inp",
                              mesh_file)

# Map the unstructured mesh with the mapping above
mesh = P4estMesh{2}(mesh_file, polydeg=5, mapping=mapping, initial_refinement_level=1)

# Refine bottom left quadrant of each tree to level 2
function refine_fn(p4est, which_tree, quadrant)
  quadrant_obj = unsafe_load(quadrant)
  if quadrant_obj.x == 0 && quadrant_obj.y == 0 && quadrant_obj.level < 3
    # return true (refine)
    return Cint(1)
  else
    # return false (don't refine)
    return Cint(0)
  end
end

# Refine recursively until each bottom left quadrant of a tree has level 2.
# The mesh will be rebalanced before the simulation starts.
refine_fn_c = @cfunction(refine_fn, Cint, (Ptr{Trixi.p4est_t}, Ptr{Trixi.p4est_topidx_t}, Ptr{Trixi.p4est_quadrant_t}))
Trixi.refine_p4est!(mesh.p4est, true, refine_fn_c, C_NULL)

semi = SemidiscretizationHyperbolic(mesh, get_time_discretization(solver),
                                    equations, initial_condition, solver,
                                    boundary_conditions=Dict(
                                      :all => BoundaryConditionDirichlet(initial_condition)
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
                                     output_directory="out_nonconforming")

amr_indicator = IndicatorLöhner(semi, variable=Trixi.density)

amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level=0,
                                      med_level=3, med_threshold=0.05,
                                      max_level=5, max_threshold=0.1)

amr_callback = AMRCallback(semi, amr_controller,
                           interval=10000000000,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

callbacks = (summary_callback, amr_callback,
             analysis_callback,
             alive_callback, save_solution)

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
