# Make sure these settings are done
# p = plot(xlims = xlimits, ylims = ylimits, aspect_ratio = :equal, showaxis = false, grid = false)

include("/home2/arpit/repositories/TrixiLW.jl/utils/plot_curved_mesh.jl")

y = sqrt(5.9^2 - 3.85^2);
plot_mesh(semi, mesh, solver, xlimits = (-2.05, 0), ylimits = (-y, y))

# Total extra elements = 246
