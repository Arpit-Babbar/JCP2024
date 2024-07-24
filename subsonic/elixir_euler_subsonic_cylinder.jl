# Channel flow around a cylinder at Mach 3
#
# Boundary conditions are supersonic Mach 3 inflow at the left portion of the domain
# and supersonic outflow at the right portion of the domain. The top and bottom of the
# channel as well as the cylinder are treated as Euler slip wall boundaries.
# This flow results in strong shock refletions / interactions as well as Kelvin-Helmholtz
# instabilities at later times as two Mach stems form above and below the cylinder.
#
# For complete details on the problem setup see Section 5.7 of the paper:
# - Jean-Luc Guermond, Murtazo Nazarov, Bojan Popov, and Ignacio Tomas (2018)
#   Second-Order Invariant Domain Preserving Approximation of the Euler Equations using Convex Limiting.
#   [DOI: 10.1137/17M1149961](https://doi.org/10.1137/17M1149961)
#
# Keywords: supersonic flow, shock capturing, AMR, unstructured curved mesh, positivity preservation, compressible Euler, 2D

using Downloads: download
using Trixi
using TrixiLW
using TrixiLW: SaveSurfacePrimitives
using LinearAlgebra

###############################################################################
# semidiscretization of the compressible Euler equations

equations = CompressibleEulerEquations2D(1.4)

@inline function initial_condition_mach01_flow(x, t, equations::CompressibleEulerEquations2D)
  # set the freestream flow parameters
  rho_freestream = 1.4
  v1 = 0.1 # This is Mach number
  v2 = 0.0
  p_freestream = 1.0

  prim = SVector(rho_freestream, v1, v2, p_freestream)
  return prim2cons(prim, equations)
end

initial_condition = initial_condition_mach01_flow

# Supersonic inflow boundary condition.
# Calculate the boundary flux entirely from the external solution state, i.e., set
# external solution state values for everything entering the domain.
@inline function boundary_condition_subsonic_constant(U_inner, f_inner, u_inner,
  outer_cache,
  normal_direction::AbstractVector, x, t, dt,
  surface_flux_function, equations::CompressibleEulerEquations2D,
  dg, time_discretization, scaling_factor = 1.0)

  u_boundary = initial_condition_mach01_flow(x, t, equations)

  return Trixi.flux_hll(u_inner, u_boundary, normal_direction, equations)
end

boundary_conditions = Dict(
  # :Bottom => TrixiLW.slip_wall_approximate,
  :Bottom => boundary_condition_subsonic_constant,
  :circle => TrixiLW.slip_wall_approximate,
  # :Top => TrixiLW.slip_wall_approximate,
  :Top => boundary_condition_subsonic_constant,
  # :Right => TrixiLW.boundary_condition_outflow,
  :Right => boundary_condition_subsonic_constant,
  :Left => boundary_condition_subsonic_constant
  )
@inline function slip_wall_cylinder(U_inner, f_inner, u_inner,
  outer_cache,
  normal_direction::AbstractVector,
  x, t, dt,
  surface_flux_function,
  equations::CompressibleEulerEquations2D,
  dg, time_discretization, scaling_factor = 1)

  return TrixiLW.slip_wall_approximate(U_inner, f_inner, u_inner, outer_cache, normal_direction,
  x, t, dt, surface_flux_function, equations, dg, time_discretization, scaling_factor)
end

boundary_conditions_channel = Dict(
  :Bottom => boundary_condition_subsonic_constant,
  :Circle => slip_wall_cylinder,
  :Top => boundary_condition_subsonic_constant,
  # :Bottom => slip_wall_approximate,
  # :Circle => slip_wall_approximate,
  # :Top => slip_wall_approximate,
  :Right => boundary_condition_subsonic_constant,
  :Left => boundary_condition_subsonic_constant)

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
  reconstruction = TrixiLW.FirstOrderReconstruction(),
  # reconstruction=TrixiLW.MUSCLReconstruction()
  # reconstruction=TrixiLW.MUSCLHancockReconstruction()
)

volume_integral = TrixiLW.VolumeIntegralFR(TrixiLW.LW())

solver = DGSEM(polydeg=polydeg, surface_flux=surface_flux,
  volume_integral=volume_integral)

# Get the unstructured quad mesh from a file (downloads the file if not available locally)
default_mesh_file = joinpath(@__DIR__, "box_circle.inp")
# default_mesh_file = joinpath(@__DIR__, "abaqus_cylinder_in_channel.inp")
# isfile(default_mesh_file) || download("https://gist.githubusercontent.com/andrewwinters5000/a08f78f6b185b63c3baeff911a63f628/raw/addac716ea0541f588b9d2bd3f92f643eb27b88f/abaqus_cylinder_in_channel.inp",
#   default_mesh_file)
mesh_file = default_mesh_file

mesh = P4estMesh{2}(mesh_file, initial_refinement_level=0)

# Running with initial_refinement_level = 4 without amr gives perfect symmetry

semi = TrixiLW.SemidiscretizationHyperbolic(mesh, get_time_discretization(solver),
  equations, initial_condition, solver, boundary_conditions=boundary_conditions)

###############################################################################
# ODE solvers

tspan = (0.0, 40.0) # Running till 40 without AMR gives perfect symmetry
lw_update = TrixiLW.semidiscretize(semi, get_time_discretization(solver), tspan);

# Callbacks

summary_callback = SummaryCallback()

analysis_interval = 8000
analysis_callback = AnalysisCallback(semi, interval=analysis_interval,
                                     output_directory = "analysis_results", save_analysis = true,
                                     analysis_integrals = (SaveSurfacePrimitives(cache -> semi.boundary_conditions.boundary_indices[1]),)
                                     )

alive_callback = AliveCallback(analysis_interval=analysis_interval)

save_solution = SaveSolutionCallback(interval=8000,
  save_initial_solution=true,
  save_final_solution=true,
  solution_variables=cons2prim)

#
cylinder_center = [0.0, 0.0]
cylinder_radius = 4.0
amr_indicator = RadialIndicator(cylinder_center, 10.0*cylinder_radius)

@inline function mach_number(u, equations)
  rho, v1, v2, p = cons2prim(u, equations)
  return sqrt(v1^2+v2^2)/sqrt(equations.gamma * p / rho)
end
# amr_indicator = IndicatorLöhner(semi, variable=mach_number)

amr_controller = ControllerThreeLevel(semi, amr_indicator,
                                      base_level=0,
                                      med_level = 1, med_threshold=0.0001, # It is a zero-one indicator
                                      max_level = 1, max_threshold=0.001)

amr_callback = AMRCallback(semi, amr_controller,
                           interval=1000000,
                           adapt_initial_condition=true,
                           adapt_initial_condition_only_refine=true)

summary_callback = SummaryCallback()
callbacks = (  analysis_callback, alive_callback, save_solution,
               amr_callback,
	       summary_callback
            )

###############################################################################
# run the simulation
time_int_tol = 1e-6
tolerances = (; abstol=time_int_tol, reltol=time_int_tol);
dt_initial = 1e-8;
cfl_number = 0.1
sol = TrixiLW.solve_lwfr(lw_update, callbacks, dt_initial, tolerances,
  # time_step_computation = TrixiLW.Adaptive(),
  time_step_computation=TrixiLW.CFLBased(cfl_number),
  # limiters=(; stage_limiter!)
);
summary_callback() # print the timer summary

