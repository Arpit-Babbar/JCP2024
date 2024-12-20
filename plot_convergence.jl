include("mpl.jl")
using Printf
using DelimitedFiles

function error_label(error_norm)
    if error_norm in ["l2", "L2"]
        return "\$L^2\$ error"
    elseif error_norm in ["l1", "L1"]
        return "\$L^1\$ error"
    end
    @assert false
end

function format_with_powers(y, _)
    if y > 1e-4
        y = Int64(y)
        return "\$$y^2\$"
    else
        return @sprintf "%.4E" y
    end
end

function set_ticks!(ax, log_sub, ticks_formatter; dim = 2)
    # Remove scientific notation and set xticks
    # https://stackoverflow.com/a/49306588/3904031

    function anonymous_formatter(y, _)
        if y > 1e-4
            # y_ = parse(Int64, y)
            y_ = Int64(y)
            if dim == 2
                return "\$$y_^2\$"
            else
                return "\$$y_\$"
            end
        else
            return @sprintf "%.4E" y
        end
    end

    formatter = plt.matplotlib.ticker.FuncFormatter(ticks_formatter)
    # (y, _) -> format"{:.4g}".format(int(y)) ) # format"{:.4g}".format(int(y)))
    # https://stackoverflow.com/questions/30887920/how-to-show-minor-tick-labels-on-log-scale-with-plt.matplotlib
    x_major = plt.matplotlib.ticker.LogLocator(base = 2.0, subs = (log_sub,),
                                               numticks = 20) # ticklabels at 2^i*log_sub
    x_minor = plt.matplotlib.ticker.LogLocator(base = 2.0,
                                               subs = LinRange(1.0, 9.0, 9) * 0.1,
                                               numticks = 10)
    #  Used to manipulate tick labels. See help(plt.matplotlib.ticker.LogLocator) for details)
    ax.xaxis.set_major_formatter(anonymous_formatter)
    ax.xaxis.set_minor_formatter(plt.matplotlib.ticker.NullFormatter())
    ax.xaxis.set_major_locator(x_major)
    ax.xaxis.set_minor_locator(x_minor)
    ax.tick_params(axis = "both", which = "major")
    ax.tick_params(axis = "both", which = "minor")
end

function add_theo_factors!(ax, ncells, error, degree, i,
    theo_factor_even, theo_factor_odd; show_label = true)
    if degree isa Int64
        d = degree
    else
        d = parse(Int64, degree)
    end
    min_y = minimum(error[1:(end-1)])
    @show error, min_y
    @show log.(error[1:end-1] ./ error[2:end]) / log(2)
    xaxis = ncells[(end-1):end]
    slope = d + 1
    if iseven(slope)
        theo_factor = theo_factor_even
    else
        theo_factor = theo_factor_odd
    end
    y0 = theo_factor * min_y
    y = (1.0 ./ xaxis) .^ slope * y0 * xaxis[1]^slope
    markers = ["", "o", "*", "^"]
    # if i == 1
    if show_label
        label = "\$ O(M^{-$(d + 1)})\$"
    else
        label = ""
    end
    ax.loglog(xaxis, y, label = label, linestyle="--",
        marker=markers[i], c="grey",
        fillstyle="none")
    # else
    # ax.loglog(xaxis,y, linestyle = "--", c = "grey")
    # end
end

function plot_python_ndofs_vs_y(files::Vector{String}, labels::Vector{String},
    degrees::Vector{Int}; y_indices = 2 * ones(Int64, length(files)),
    saveto,
    theo_factor_even=0.8, theo_factor_odd=0.8,
    title=nothing, log_sub="2.5",
    error_norm="l2",
    ticks_formatter=format_with_powers,
    figsize=(6.4, 4.8))
    # @assert error_type in ["l2","L2"] "Only L2 error for now"
    fig_error, ax_error = plt.subplots(figsize=figsize)
    colors = ["orange", "royalblue", "green", "m", "c", "y", "k"]
    markers = ["D", "o", "*", "^"]
    @assert length(files) == length(labels)
    n_plots = length(files)

    for i in 1:n_plots
        degree = degrees[i]
        data = readdlm(files[i])
        marker = markers[i]
        # TODO - This sqrt shows really bad foresight!
        ax_error.loglog(data[:, 1], data[:, y_indices[i]], marker=marker, c=colors[i],
            mec="k", fillstyle="none", label=labels[i])
    end

    done_degrees = Vector{Int64}()
    for i in eachindex(degrees) # Assume degrees are not repeated
        data = readdlm(files[i])
        degree = degrees[i]
        if degree in done_degrees
            show_label = false
        else
            show_label = true
        end
        push!(done_degrees, degree)
        add_theo_factors!(ax_error, data[:, 1], data[:, y_indices[i]], degree, i,
            theo_factor_even, theo_factor_odd, show_label = show_label)
    end
    ax_error.set_xlabel("Number of elements")
    ax_error.set_ylabel(error_label(error_norm))

    set_ticks!(ax_error, log_sub, ticks_formatter; dim=2)

    ax_error.grid(true, linestyle="--")

    if title !== nothing
        ax_error.set_title(title)
    end
    ax_error.legend()

    fig_error.savefig("$saveto.pdf")
    # fig_error.savefig("$saveto.png")

    return fig_error
end

files = ["isentropic/errors_N3.txt", "isentropic/errors_N3.txt", "isentropic/errors_N3.txt", "isentropic/errors_N3.txt"]
y_indices = [2, 3, 4, 5]
labels = ["\$ \\rho \$", "\$ \\rho v_1 \$", "\$ \\rho v_2 \$", "\$ E \$"]
degrees = [3, 3, 3, 3]

plot_python_ndofs_vs_y(files, labels, degrees, saveto = "convergence", log_sub = "2.0",
                       figsize = (6.0, 6.5), title = "Degree \$ N = 3 \$", y_indices = y_indices)
