using ROOTPlots

include(joinpath(@__DIR__, "..", "examples", "synthetic_control.jl"))

"""Save a synthetic three-panel control figure using native ROOT only."""
function main(args)
    length(args) <= 1 || error("usage: demo.jl [output-stem]")
    root_available() || error("native ROOT is required; install ROOT and check ROOT_CONFIG or PATH")
    outstem = isempty(args) ? joinpath(@__DIR__, "..", "examples", "output", "control_demo") : args[1]
    cmps, specs = demo_panels()
    result = render_root(cmps, specs, outstem; format=:both)
    println("Backend: native ROOT")
    println("PNG: ", result.png)
    println("PDF: ", result.pdf)
end

main(ARGS)
