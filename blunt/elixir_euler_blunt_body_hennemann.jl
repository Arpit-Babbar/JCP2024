using Trixi, TrixiLW
using Plots
plotlyjs()
using Trixi: transfinite_mapping

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

function initial_value_blunt_body(x, t, equations)
   return prim2cons((1.4, 4.0, 0.0, 1.0), equations)
end
polydeg = 4
basis = LobattoLegendreBasis(polydeg)
initial_condition = initial_value_blunt_body
surface_flux=flux_lax_friedrichs
shock_indicator = IndicatorHennemannGassner(equations, basis,
					   alpha_max = 1.0,
					   alpha_min = 0.001,
					   alpha_smooth = true,
					   variable = density_pressure)
volume_integral = TrixiLW.VolumeIntegralFRShockCapturing(shock_indicator;
volume_flux_fv=surface_flux, reconstruction=TrixiLW.FirstOrderReconstruction())

solver = DGSEM(polydeg=polydeg,surface_flux = surface_flux ,
               volume_integral=volume_integral)

@inline function boundary_condition_supersonic_inflow(U_inner, f_inner, u_inner,
   outer_cache, normal_direction::AbstractVector, x, t, dt,
   surface_flux_function, equations::CompressibleEulerEquations2D,
   dg, time_discretization, scaling_factor = 1)
   u_boundary = initial_value_blunt_body(x, t, equations)
   flux = Trixi.flux(u_boundary, normal_direction, equations)
   return flux
end

@inline function boundary_condition_subsonic_constant(U_inner, f_inner, u_inner,
   outer_cache,
   normal_direction::AbstractVector, x, t, dt,
   surface_flux_function, equations::CompressibleEulerEquations2D,
   dg, time_discretization, scaling_factor = 1.0)

   u_boundary = initial_value_blunt_body(x, t, equations)

   return Trixi.flux_hll(u_inner, u_boundary, normal_direction, equations)
end

boundary_conditions = Dict(:body=>TrixiLW.slip_wall_approximate,
                       :outer_arc => boundary_condition_subsonic_constant
                       )


mesh_file = "blunt_body.inp"
mesh = P4estMesh{2}(mesh_file, initial_refinement_level=0)

cfl_number = 0.05
semi = TrixiLW.SemidiscretizationHyperbolic(mesh, get_time_discretization(solver),
 equations, initial_condition, solver,
 boundary_conditions = boundary_conditions
 )

stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds=(5.0e-7, 1.0e-6),
 variables=(Trixi.pressure, Trixi.density))

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 10.0) # Use tfinal = 1.0, with cfl_number = 0.05 for benchmarking
# tspan = (0.0, 0.0) # Use tfinal = 1.0, with cfl_number = 0.05 for benchmarking
lw_update = TrixiLW.semidiscretize(semi, get_time_discretization(solver), tspan);

analysis_interval = 500
analysis_callback = AnalysisCallback(semi, interval=analysis_interval)

alive_callback = AliveCallback(analysis_interval=analysis_interval)

visualization_callback = VisualizationCallback(interval=100,
   save_initial_solution=true,
   save_final_solution=true)

save_restart = SaveRestartCallback(interval=10000,
                                   save_final_restart=true)
save_solution = SaveSolutionCallback(interval=500,
                                     save_initial_solution=true,
                                     save_final_solution=true,
                                     solution_variables=cons2prim)

amr_indicator = IndicatorLöhner(semi, variable=Trixi.density)

amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level=0,
                                      med_level=1, med_threshold=0.05,
                                      # med_level=3, med_threshold=0.05,
                                      max_level=2, max_threshold=0.1
                                      # max_level=4, max_threshold=0.13
				      )

amr_callback = AMRCallback(semi, amr_controller,
                           interval=1,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

summary_callback = SummaryCallback()
callbacks = ( analysis_callback, alive_callback,
               save_solution,
               amr_callback, summary_callback
            );

###############################################################################
# run the simulation

time_int_tol = 1e-6
tolerances = (;abstol = time_int_tol, reltol = time_int_tol);
dt_initial = 1e-8;
sol = TrixiLW.solve_lwfr(lw_update, callbacks, dt_initial, tolerances,
                     time_step_computation = TrixiLW.Adaptive(),
                     # time_step_computation = TrixiLW.CFLBased(cfl_number),
                      limiters = (;stage_limiter!)
                      );
summary_callback() # print the timer summary

