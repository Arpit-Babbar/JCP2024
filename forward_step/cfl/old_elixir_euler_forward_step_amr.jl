using Downloads: download
using TrixiLW
using Trixi
using LinearAlgebra

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

"""
    initial_condition_mach3_flow(x, t, equations::CompressibleEulerEquations2D)

Compressible Euler setup for a Mach 3 wind tunnel flow with a forward facing step.
Results in strong shock interactions as well as Kelvin-Helmholtz instabilities at later times.
See Section IV b on the paper below for details.

- Paul Woodward and Phillip Colella (1984)
  The Numerical Simulation of Two-Dimensional Fluid Flows with Strong Shocks.
  [DOI: 10.1016/0021-9991(84)90142-6](https://doi.org/10.1016/0021-9991(84)90142-6)
"""
@inline function initial_condition_mach3_flow(x, t, equations::CompressibleEulerEquations2D)
   # set the freestream flow parameters
   rho_freestream = 1.4
   v1 = 3.0
   v2 = 0.0
   p_freestream = 1.0

   prim = SVector(rho_freestream, v1, v2, p_freestream)
   return prim2cons(prim, equations)
end

initial_condition = initial_condition_mach3_flow

@inline function boundary_condition_inflow_forward_step(U_inner, f_inner, u_inner,
   outer_cache,
   normal_direction::AbstractVector, x, t, dt,
   surface_flux_function, equations::CompressibleEulerEquations2D,
   dg, time_discretization, scaling_factor = 1)
   u_outer = initial_condition_mach3_flow(x, t, equations)
   flux = Trixi.flux(u_outer, normal_direction, equations)

   return flux
end

@inline function boundary_condition_outflow(U_inner, f_inner, u_inner,
   outer_cache,
   normal_direction::AbstractVector, x, t, dt,
   surface_flux_function, equations::CompressibleEulerEquations2D,
   dg, time_discretization, scaling_factor = 1)
   # flux = Trixi.flux(u_inner, normal_direction, equations)

   # return flux
   return f_inner
end

@inline function slip_wall_approximate(U_inner, f_inner, u_inner,
   outer_cache,
   normal_direction::AbstractVector,
   x, t, dt,
   surface_flux_function,
   equations::CompressibleEulerEquations2D,
   dg, time_discretization, scaling_factor = 1)

   # TRIXI WAY!
   # return boundary_condition_slip_wall(u_inner, normal_direction, x, t, surface_flux_function, equations)
   norm_ = norm(normal_direction)
   # Normalize the vector without using `normalize` since we need to multiply by the `norm_` later
   normal = normal_direction / norm_

   # rotate the internal solution state
   u_local = Trixi.rotate_to_x(u_inner, normal, equations)

   # compute the primitive variables
   rho_local, v_normal, v_tangent, p_local = cons2prim(u_local, equations)

   # Get the solution of the pressure Riemann problem
   # See Section 6.3.3 of
   # Eleuterio F. Toro (2009)
   # Riemann Solvers and Numerical Methods for Fluid Dynamics: A Pratical Introduction
   # [DOI: 10.1007/b79761](https://doi.org/10.1007/b79761)
   if v_normal <= 0.0
      sound_speed = sqrt(equations.gamma * p_local / rho_local) # local sound speed
      a_ = 1.0 + 0.5 * (equations.gamma - 1) * v_normal / sound_speed
      if a_ >= 0.0
         p_star = p_local * (a_)^(2.0 * equations.gamma * equations.inv_gamma_minus_one)
      else
         u_outer = SSFR.get_reflection(u_inner, normal_direction, equations)

         flux = surface_flux_function(u_inner, u_outer, normal_direction, equations)

         return flux
      end
   else # v_normal > 0.0
      # @show v_normal, p_local
      A = 2.0 / ((equations.gamma + 1) * rho_local)
      B = p_local * (equations.gamma - 1) / (equations.gamma + 1)
      p_star = p_local + 0.5 * v_normal / A * (v_normal + sqrt(v_normal^2 + 4.0 * A * (p_local + B)))
   end

   # For the slip wall we directly set the flux as the normal velocity is zero
   return SVector(zero(eltype(u_inner)),
      p_star * normal[1],
      p_star * normal[2],
      zero(eltype(u_inner))) * norm_
end

boundary_conditions = Dict(:Bottom => slip_wall_approximate,
   :Step_Front => slip_wall_approximate,
   :Step_Top => slip_wall_approximate,
   :Top => slip_wall_approximate,
   :Right => boundary_condition_outflow,
   :Left => boundary_condition_inflow_forward_step)

# boundary_conditions = Dict(:Bottom => TrixiLW.boundary_condition_slip_wall_horizontal,
#    :Step_Front => TrixiLW.boundary_condition_slip_wall_vertical,
#    :Step_Top => TrixiLW.boundary_condition_slip_wall_horizontal,
#    :Top => TrixiLW.boundary_condition_slip_wall_horizontal,
#    :Right => TrixiLW.boundary_condition_outflow,
#    :Left => boundary_condition_inflow_forward_step)

surface_flux = flux_lax_friedrichs

polydeg = 4
basis = LobattoLegendreBasis(polydeg)
shock_indicator = IndicatorHennemannGassner(equations, basis,
   alpha_max=1.0,
   alpha_min=0.001,
   alpha_smooth=true,
   variable=density_pressure)
volume_integral = TrixiLW.VolumeIntegralFRShockCapturing(
   shock_indicator;
   volume_flux_fv=surface_flux,
   reconstruction=TrixiLW.FirstOrderReconstruction()
   # reconstruction=TrixiLW.MUSCLReconstruction()
   # reconstruction=TrixiLW.MUSCLHancockReconstruction()
)

solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux,
   volume_integral=volume_integral)

# Get the unstructured quad mesh from a file (downloads the file if not available locally)
default_mesh_file = joinpath(@__DIR__, "abaqus_forward_step.inp")
isfile(default_mesh_file) || download("https://gist.githubusercontent.com/andrewwinters5000/b346ee6aa5446687f128eab8b37d52a7/raw/cd1e1d43bebd8d2631a07caec45585ec8456ca4c/abaqus_forward_step.inp",
   default_mesh_file)
mesh_file = default_mesh_file

mesh = P4estMesh{2}(mesh_file)

semi = TrixiLW.SemidiscretizationHyperbolic(mesh, get_time_discretization(solver),
  equations, initial_condition, solver, boundary_conditions=boundary_conditions)

###############################################################################
# ODE solvers, callbacks etc.

tspan = (0.0, 3.0)
lw_update = TrixiLW.semidiscretize(semi, get_time_discretization(solver), tspan);

summary_callback = SummaryCallback()

analysis_interval = 10000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
   extra_analysis_integrals=(entropy,))

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=20000,
   save_initial_solution=true,
   save_final_solution=true,
   solution_variables=cons2prim)

amr_indicator = IndicatorLöhner(semi, variable=Trixi.density)

amr_controller = ControllerThreeLevel(semi, amr_indicator,
   base_level=0,
   med_level=2, med_threshold=0.05,
   max_level=5, max_threshold=0.1)

amr_callback = AMRCallback(semi, amr_controller,
   interval=100,
   adapt_initial_condition=true,
   adapt_initial_condition_only_refine=true)

callbacks = (analysis_callback, alive_callback, save_solution, amr_callback)

# positivity limiter necessary for this example with strong shocks
stage_limiter! = PositivityPreservingLimiterZhangShu(thresholds=(5.0e-6, 5.0e-6),
   variables=(Trixi.density, pressure))

###############################################################################
# run the simulation
time_int_tol = 1e-5
tolerances = (; abstol=time_int_tol, reltol=time_int_tol);
dt_initial = 1e-6;
cfl_number = 0.18 # crashes for 0.19
sol = TrixiLW.solve_lwfr(lw_update, callbacks, dt_initial, tolerances,
  # time_step_computation = TrixiLW.Adaptive(),
  time_step_computation=TrixiLW.CFLBased(cfl_number),
  limiters=(; stage_limiter!)
);
summary_callback() # print the timer summary
