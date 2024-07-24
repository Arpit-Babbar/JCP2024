Running this code requires julia. I recommend using `juliaup` to install it, see https://julialang.org/downloads/. Once you have julia installed, enter the current directory in your terminal and enter the julia REPL by using the command `julia`.

```julia
julia> import Pkg
julia> Pkg.add("https://github.com/Arpit-Babbar/TrixiLW.jl")
```
The above command requires access to `TrixiLW.jl`, which is currently private. Kindly email me to get access to it. For now, this repository shows the elixir files. Now, you are ready to run the any code, e.g., by entering
```julia
julia> include("forward_step/error/elixir_euler_forward_step_amr.jl")
```
