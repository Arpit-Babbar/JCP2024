using Trixi: trixi_include
using DelimitedFiles

n_levels = 6
min_level = 0
max_level = n_levels - 1

levels = min_level:max_level

trees_per_dimension = (8,8)

l2_errors = zeros(n_levels, 5)
nxs = zeros(n_levels)

println("Running isentropic convergence")
for degree in 3:3
  for (i, level) in enumerate(levels)
    trixi_include(joinpath(@__DIR__,"elixir_euler_isentropic_hennemann.jl"),
              initial_refinement_level = level,
      cfl_number = 0.25,
      polydeg = degree)
    for n in 1:4
      l2_errors[i,n] = analysis_callback(sol).l2[n]
    end
    nxs[i] = 16 * 2^level
  end

  data = hcat(nxs, l2_errors)
  writedlm(joinpath(@__DIR__,"errors_N$degree.txt"), data)
end
